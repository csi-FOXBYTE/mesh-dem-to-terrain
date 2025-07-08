import { Document, Logger } from "@gltf-transform/core";
import { dedup, flatten, join, reorder, weld } from "@gltf-transform/functions";
import Delaunator from "delaunator";
import { readFile } from "fs/promises";
import _ from "lodash";
import { MeshoptEncoder } from "meshoptimizer";
import RBush from "rbush";
import { Box3, Vector3 } from "three";
import { parentPort } from "worker_threads";
import {
  betterWeld,
  cut,
  expandWithExternalPolygons,
  flipWindingOrder,
  getSkirtIndicesFlat,
  resampleGrid,
  splitSharpEdges,
} from "./functions.js";
import { QuantizedMesh } from "./quantizedMesh.js";
import { Tile } from "./types.js";
import { WorkerPayloads } from "./workerTypes.js";
import path from "path";

Logger.DEFAULT_INSTANCE = new Logger(Logger.Verbosity.SILENT);

if (!parentPort) throw new Error("Is not being called in a worker context!");

let tree: RBush<Tile> | null = null;
let outputFolder: string | null = null;
let skipLevel: number = 10;

parentPort.on("message", async (message: WorkerPayloads) => {
  try {
    switch (message.type) {
      case "init":
        tree = new RBush<Tile>().fromJSON(message.tree);
        parentPort!.postMessage("ok");
        outputFolder = message.outputFolder;
        skipLevel = message.skipLevel;
        break;
      case "work":
        if (!tree) throw new Error("Worker is not initialized!");
        if (!outputFolder) throw new Error("No outputfolder set!");

        const [mx0, my0, mx1, my1] = tmsTileToEPSG3857(
          message.x,
          message.y,
          message.zoom
        );

        const document = new Document();

        const buffer = document.createBuffer();
        const root = document.getRoot();
        const scene = document.createScene();
        root.setDefaultScene(scene);

        const material = document.createMaterial();

        const dx = (mx1 - mx0) * 0.1;
        const dy = (my1 - my0) * 0.1;

        const tileBBox = new Box3(
          new Vector3(mx0, my0, -2000),
          new Vector3(mx1, my1, 2000)
        );

        const searched =
          message.zoom < skipLevel
            ? []
            : tree.search({
                minX: mx0 - dx,
                minY: my0 - dy,
                maxX: mx1 + dx,
                maxY: my1 + dy,
              });

        for (const { path } of searched) {
          const { gt, indices, vertices } = await deserializeIndicesAndVertices(
            path
          );

          const primitive = document.createPrimitive();

          const positionAccessor = document
            .createAccessor()
            .setBuffer(buffer)
            .setArray(vertices)
            .setType("VEC3");

          const indicesAccessor = document
            .createAccessor()
            .setBuffer(buffer)
            .setArray(indices)
            .setType("SCALAR");

          primitive.setAttribute("POSITION", positionAccessor);
          primitive.setIndices(indicesAccessor);

          primitive.setMaterial(material);

          const node = document.createNode();
          scene.addChild(node);

          node.setTranslation(gt.toArray());

          const mesh = document.createMesh();
          node.setMesh(mesh);

          mesh.addPrimitive(primitive);
        }

        if (searched.length !== 0) {
          const gt = new Vector3(
            ...(root.listNodes()?.[0]?.getTranslation() ?? [])
          );

          if (message.zoom < 16) {
            await document.transform(
              flipWindingOrder(),
              dedup(),
              flatten(),
              join(),
              weld({}),
              resampleGrid({
                minX: mx0 - gt.x,
                minY: my0 - gt.y,
                maxX: mx1 - gt.x,
                maxY: my1 - gt.y,
                offset: -100,
              }),
              splitSharpEdges(),
              reorder({
                encoder: MeshoptEncoder,
              })
            );
          } else {
            await document.transform(
              flipWindingOrder(),
              dedup(),
              flatten(),
              join(),
              betterWeld(),
              expandWithExternalPolygons({
                minX: mx0 - gt.x,
                minY: my0 - gt.y,
                maxX: mx1 - gt.x,
                maxY: my1 - gt.y,
                offset: -100,
              }),
              cut({
                minX: mx0 - gt.x,
                minY: my0 - gt.y,
                maxX: mx1 - gt.x,
                maxY: my1 - gt.y,
              }),
              weld({}),
              splitSharpEdges(),
              reorder({
                encoder: MeshoptEncoder,
              })
            );
          }
        }

        const vertexCounts = root
          .listMeshes()
          .flatMap((mesh) => mesh.listPrimitives())
          .flatMap(
            (prim) => prim.getAttribute("POSITION")?.getArray()?.byteLength
          );

        const sum = _.sum(vertexCounts);

        let vertices: number[] = [];
        let indices: number[] = [];

        if (sum === 0) {
          const stepX = (mx1 - mx0) / 31;
          const stepY = (my1 - my0) / 31;

          for (let y = 0; y < 32; y++) {
            for (let x = 0; x < 32; x++) {
              vertices.push(mx0 + stepX * x, my0 + stepY * y, -100);
            }
          }

          const { triangles } = new Delaunator(
            vertices.filter((_, i) => (i + 1) % 3 !== 0)
          );

          triangles.reverse();

          const [remap, unique] = MeshoptEncoder.reorderMesh(
            triangles,
            true,
            false
          );

          const newVertices: number[] = new Array(unique * 3);
          const remapMap = new Map<number, number>();

          for (let i = 0; i < vertices.length / 3; i++) {
            if (remap[i] != 0xffffffff) {
              const pIndex = i * 3;
              const rIndex = remap[i] * 3;
              newVertices[rIndex] = vertices[pIndex];
              newVertices[rIndex + 1] = vertices[pIndex + 1];
              newVertices[rIndex + 2] = vertices[pIndex + 2];

              remapMap.set(i, remap[i]);
            }
          }

          vertices = newVertices;

          indices = Array.from(triangles);
        } else {
          vertices = Array.from(
            root
              .listMeshes()
              .flatMap((mesh) => mesh.listPrimitives())[0]
              .getAttribute("POSITION")!
              .getArray()!
          );

          const gt = new Vector3(
            ...(root.listNodes()?.[0]?.getTranslation() ?? [])
          );

          for (let i = 0; i < vertices.length; i += 3) {
            vertices[i + 0] += gt.x;
            vertices[i + 1] += gt.y;
            vertices[i + 2] += gt.z;
          }

          indices = Array.from(
            root
              .listMeshes()
              .flatMap((mesh) => mesh.listPrimitives())[0]
              .getIndices()!
              .getArray()!
          );
        }

        const { east, north, south, west } = getSkirtIndicesFlat(
          vertices,
          indices,
          3,
          1e-6
        );

        const quantizedMesh = await new QuantizedMesh(
          vertices,
          indices,
          east,
          west,
          south,
          north,
          tileBBox
        ).serialize();

        parentPort!.postMessage({
          path: path.join(
            outputFolder,
            message.zoom.toString(),
            message.x.toString(),
            `${message.y}.terrain`
          ),
          file: quantizedMesh,
        });
    }
  } catch (e) {
    console.error(e);
  }
});

async function deserializeIndicesAndVertices(input: string) {
  const file = await readFile(input);

  const gt = new Vector3();

  const x = file.readDoubleLE(0);
  const y = file.readDoubleLE(8);
  const z = file.readDoubleLE(16);

  const indexStart = 24;

  const indicesLength = file.readInt32LE(indexStart);

  const indices = new Uint32Array(indicesLength);

  const isBigIndex = indicesLength > 2 ** 16;
  const indexStep = isBigIndex ? 4 : 2;

  if (isBigIndex) {
    for (let i = 0; i < indicesLength; i++) {
      indices[i] = file.readUint32LE(indexStart + 4 + i * indexStep);
    }
  } else {
    for (let i = 0; i < indicesLength; i++) {
      indices[i] = file.readUint16LE(indexStart + 4 + i * indexStep);
    }
  }

  const vertexStart = indexStart + 4 + indicesLength * indexStep;

  const verticesLength = file.readInt32LE(vertexStart);

  const vertices = new Float32Array(verticesLength);

  for (let i = 0; i < verticesLength; i++) {
    vertices[i] = file.readFloatLE(vertexStart + 4 + i * 4);
  }

  gt.set(x, y, z);

  return { gt, indices, vertices };
}

/**
 * Convert a TMS tile coordinate (x, y, z) to EPSG:3857 bounds in metres.
 *
 * TMS uses origin at bottom‐left, tileSize pixels per tile (default 256).
 * EPSG:3857 (Web Mercator) spans [−originShift, +originShift] in metres.
 */
export function tmsTileToEPSG3857(
  x: number,
  y: number,
  z: number,
  tileSize: number = 256
) {
  const RADIUS = 6378137; // earth radius in metres
  const ORIGIN_SHIFT = (2 * Math.PI * RADIUS) / 2; // half the world span
  const INITIAL_RESOLUTION = (2 * Math.PI * RADIUS) / tileSize;
  const resolution = INITIAL_RESOLUTION / Math.pow(2, z);

  const minX = x * tileSize * resolution - ORIGIN_SHIFT;
  const minY = y * tileSize * resolution - ORIGIN_SHIFT;
  const maxX = (x + 1) * tileSize * resolution - ORIGIN_SHIFT;
  const maxY = (y + 1) * tileSize * resolution - ORIGIN_SHIFT;

  return [minX, minY, maxX, maxY] as const;
}
