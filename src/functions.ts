import { Transform, Document, Node, Primitive } from "@gltf-transform/core";
import { polygonclip } from "./lineclip/index.js";
import { triangulate3DPolygon } from "./helpers.js";
import { compactPrimitive } from "@gltf-transform/functions";
import Delaunator from "delaunator";
import { Box3, Vector3 } from "three";
import FlatBush from "flatbush";
import earcut from "earcut";
import _ from "lodash";

export function flipWindingOrder(): Transform {
  return (document: Document) => {
    const root = document.getRoot();

    const primitives = root.listMeshes().flatMap((m) => m.listPrimitives());

    for (const primitive of primitives) {
      const accessor = primitive.getIndices();
      const idx = accessor?.getArray();
      if (!idx) return;

      for (let i = 0; i < idx.length; i += 3) {
        // [i0, i1, i2] → [i0, i2, i1]
        const tmp = idx[i + 1];
        idx[i + 1] = idx[i + 2];
        idx[i + 2] = tmp;
      }

      accessor?.setArray(idx);
    }
  };
}

export function resampleGrid({
  maxX,
  maxY,
  minX,
  minY,
  offset,
}: {
  minX: number;
  maxX: number;
  minY: number;
  maxY: number;
  offset: number;
}): Transform {
  return (document) => {
    for (const [node, primitive] of document
      .getRoot()
      .listNodes()
      .flatMap((node) => {
        return (node
          .getMesh()
          ?.listPrimitives()
          .map((primitive) => [node, primitive]) ?? []) as [Node, Primitive][];
      })) {
      const indices = primitive.getIndices()?.getArray();
      const positions = primitive.getAttribute("POSITION")?.getArray();

      if (!indices || !positions) continue;

      const vertices: [number, number, number][] = new Array(65 * 65)
        .fill(0)
        .map(() => [0, 0, 0]); // left to right, up to down

      const zs = new Float32Array(indices.length / 3);

      const tree = new FlatBush(indices.length / 3);

      const v0 = new Vector3();
      const v1 = new Vector3();
      const v2 = new Vector3();

      const triBox = new Box3();

      let t = 0;
      for (let i = 0; i < indices.length; i += 3, t++) {
        const i0 = indices[i] * 3;
        const i1 = indices[i + 1] * 3;
        const i2 = indices[i + 2] * 3;

        v0.set(positions[i0], positions[i0 + 1], positions[i0 + 2]);
        v1.set(positions[i1], positions[i1 + 1], positions[i1 + 2]);
        v2.set(positions[i2], positions[i2 + 1], positions[i2 + 2]);

        triBox.makeEmpty();

        triBox.setFromPoints([v0, v1, v2]);

        tree.add(triBox.min.x, triBox.min.y, triBox.max.x, triBox.max.y);

        zs[t] = (v0.z + v1.z + v2.z) / 3;
      }

      tree.finish();

      const sizeX = maxX - minX;
      const sizeY = maxY - minY;

      for (let y = 0; y < 65; y++) {
        for (let x = 0; x < 65; x++) {
          const index = x + y * 65;

          const calculatedMinX = (x / 64) * sizeX + minX;
          const calculatedMinY = (y / 64) * sizeY + minY;

          const calculatedMaxX = ((x + 1) / 64) * sizeX + minX;
          const calculatedMaxY = ((y + 1) / 64) * sizeY + minY;

          let finalZ = 0;

          const foundEntries = tree.search(
            calculatedMinX,
            calculatedMinY,
            calculatedMaxX,
            calculatedMaxY
          );

          if (foundEntries.length === 0) {
            vertices[index] = [
              calculatedMinX,
              calculatedMinY,
              offset - node.getTranslation()[2],
            ];
            continue;
          }

          foundEntries.forEach((entry) => (finalZ += zs[entry]));

          vertices[index] = [
            calculatedMinX,
            calculatedMinY,
            finalZ / foundEntries.length,
          ];
        }
      }

      const { triangles } = Delaunator.from(
        vertices,
        (p) => p[0],
        (p) => p[1]
      );

      for (let i = 0; i < triangles.length; i += 3) {
        // [i0, i1, i2] → [i0, i2, i1]
        const tmp = triangles[i + 1];
        triangles[i + 1] = triangles[i + 2];
        triangles[i + 2] = tmp;
      }

      primitive
        .getAttribute("POSITION")
        ?.setArray(new Float32Array(vertices.flat()));
      primitive.getIndices()?.setArray(new Uint32Array(triangles));
    }
  };
}

export function cut({
  maxX,
  maxY,
  minX,
  minY,
}: {
  minX: number;
  maxX: number;
  minY: number;
  maxY: number;
}): Transform {
  return (document) => {
    const root = document.getRoot();

    const primitives = root.listMeshes().flatMap((m) => m.listPrimitives());

    const bbox = [minX, minY, maxX, maxY] as [number, number, number, number];

    for (const primitive of primitives) {
      const vertices = primitive.getAttribute("POSITION")?.getArray();
      const indices = primitive.getIndices()?.getArray();

      if (!vertices || !indices) continue;

      const discardedTriangleIndicesSet = new Set<number>();

      const createdVertices: number[] = [];
      const createdIndices: number[] = [];

      for (let i = 0; i < indices.length; i += 3) {
        const i0 = indices[i] * 3;
        const i1 = indices[i + 1] * 3;
        const i2 = indices[i + 2] * 3;

        const originalPolygon = [
          [vertices[i0], vertices[i0 + 1], vertices[i0 + 2]],
          [vertices[i1], vertices[i1 + 1], vertices[i1 + 2]],
          [vertices[i2], vertices[i2 + 1], vertices[i2 + 2]],
        ] as [number, number, number][];

        const newPolygon = polygonclip(originalPolygon, bbox);

        if (newPolygon.length === 0) {
          discardedTriangleIndicesSet.add(i);
        } else {
          discardedTriangleIndicesSet.add(i);
          const indices = triangulate3DPolygon([newPolygon]);

          if (!indices) {
            continue;
          }
          createdIndices.push(
            ...indices.map((index) => index + createdVertices.length / 3)
          );
          createdVertices.push(...newPolygon.flat());
        }
      }

      const newIndices: number[] = [];
      const newVertices: number[] = Array.from(vertices);

      for (let i = 0; i < indices.length; i += 3) {
        if (discardedTriangleIndicesSet.has(i)) continue;

        const i0 = indices[i];
        const i1 = indices[i + 1];
        const i2 = indices[i + 2];

        newIndices.push(i0, i1, i2);
      }

      primitive
        .getAttribute("POSITION")
        ?.setArray(new Float32Array(newVertices.concat(createdVertices)));
      primitive
        .getIndices()
        ?.setArray(
          new Uint32Array(
            newIndices.concat(
              createdIndices.map((index) => index + vertices.length / 3)
            )
          )
        );
      compactPrimitive(primitive);
    }
  };
}

export function expandWithExternalPolygons({
  maxX,
  maxY,
  minX,
  minY,
  offset,
}: {
  minX: number;
  maxX: number;
  minY: number;
  maxY: number;
  offset: number;
}): Transform {
  return (document) => {
    for (const [node, primitive] of document
      .getRoot()
      .listNodes()
      .flatMap((node) => {
        return (node
          .getMesh()
          ?.listPrimitives()
          .map((primitive) => [node, primitive]) ?? []) as [Node, Primitive][];
      })) {
      const posAttr = primitive.getAttribute("POSITION");
      const idxAttr = primitive.getIndices();
      const vertices = posAttr?.getArray();
      const indices = idxAttr?.getArray();

      if (!vertices || !indices) {
        continue;
      }

      // --- extract 2D rings from original mesh ---
      const verts2d: [number, number][] = [];
      for (let i = 0; i < vertices.length; i += 3) {
        verts2d.push([vertices[i], vertices[i + 1]]);
      }
      const boundaryRings: [number, number][][] = extractOuterBoundaryRings(
        verts2d,
        Array.from(indices)
      ).map((ring) => ring.map((p) => [p[0], p[1]]));

      if (isBoxInsidePolygon(boundaryRings.flat(), [minX, minY, maxX, maxY])) {
        continue;
      }

      // --- compute overall bounds and create big outer rectangle ---
      const min: [number, number] = [0, 0];
      posAttr!.getMin(min);
      const max: [number, number] = [0, 0];
      posAttr!.getMax(max);

      const tMinX = Math.min(minX, min[0]);
      const tMinY = Math.min(minY, min[1]);
      const tMaxX = Math.max(maxX, max[0]);
      const tMaxY = Math.max(maxY, max[1]);

      const outerRing: [number, number][] = [
        [tMinX * 2, tMinY * 2],
        [tMaxX * 2, tMinY * 2],
        [tMaxX * 2, tMaxY * 2],
        [tMinX * 2, tMaxY * 2],
      ];

      // --- remember how many verts we had ---
      const origVertCount = vertices.length / 3;

      // --- compute bottom face triangulation with earcut ---
      let holeStart = outerRing.length;
      const holeIndices = boundaryRings.map((r) => {
        const start = holeStart;
        holeStart += r.length;
        return start;
      });
      const all2D = outerRing.flat().concat(boundaryRings.flat().flat());
      const bottomTris = earcut(all2D, holeIndices, 2);

      // --- append all the new vertices at z = offset ---
      const bottomOuter = outerRing.flatMap((p) => [
        p[0],
        p[1],
        offset - node.getTranslation()[2],
      ]);
      const bottomHoles = boundaryRings
        .flat()
        .flatMap((p) => [p[0], p[1], offset - node.getTranslation()[2]]);
      posAttr!.setArray(
        new Float32Array([...vertices, ...bottomOuter, ...bottomHoles])
      );

      // --- build a lookup of original vertex‐indices for each ring ---
      const boundaryOrigIdx = boundaryRings.map((ring) =>
        ring.map((pt) =>
          verts2d.findIndex((v) => v[0] === pt[0] && v[1] === pt[1])
        )
      );

      // --- compute signed area of each ring (shoelace) to detect winding ---
      const ringAreas = boundaryRings.map((ring) => {
        let sum = 0;
        for (let i = 0; i < ring.length; i++) {
          const [x1, y1] = ring[i];
          const [x2, y2] = ring[(i + 1) % ring.length];
          sum += x1 * y2 - x2 * y1;
        }
        return sum / 2;
      });

      // --- emit side‐wall triangles with consistent winding ---
      const sideIndices: number[] = [];
      boundaryOrigIdx.forEach((origIdxRing, ringIdx) => {
        const isCCW = ringAreas[ringIdx] > 0;
        // bottom‐ring base index for this hole
        const offsetBefore = boundaryRings
          .slice(0, ringIdx)
          .reduce((sum, r) => sum + r.length, 0);
        const bottomStart = origVertCount + outerRing.length + offsetBefore;

        for (let i = 0; i < origIdxRing.length; i++) {
          const j = (i + 1) % origIdxRing.length;
          const oA = origIdxRing[i],
            oB = origIdxRing[j],
            bA = bottomStart + i,
            bB = bottomStart + j;

          if (!isCCW) {
            // same as before: (oA→oB→bB) and (bB→bA→oA)
            sideIndices.push(oA, oB, bB, bB, bA, oA);
          } else {
            // reversed winding: swap last two verts in each triangle
            sideIndices.push(oA, bB, oB, bB, oA, bA);
          }
        }
      });

      // --- finally, rebuild the full index buffer ---
      idxAttr!.setArray(
        new Uint32Array([
          ...indices,
          ...bottomTris.map((i) => i + origVertCount),
          ...sideIndices,
        ])
      );

      compactPrimitive(primitive);
    }
  };
}

type Point = [number, number];
type BBox = [minX: number, minY: number, maxX: number, maxY: number];

/**
 * Returns true if the given axis-aligned bounding box
 * is entirely inside the polygon.
 */
function isBoxInsidePolygon(polygon: Point[], box: BBox): boolean {
  const [minX, minY, maxX, maxY] = box;
  // the 4 corners
  const corners: Point[] = [
    [minX, minY],
    [maxX, minY],
    [maxX, maxY],
    [minX, maxY],
  ];

  // ray-cast test for one point
  const inside = (pt: Point) => {
    let inside = false;
    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
      const [xi, yi] = polygon[i];
      const [xj, yj] = polygon[j];
      const [x, y] = pt;
      // does edge (j→i) straddle the horizontal ray? if so, does it cross to the right of pt?
      if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) {
        inside = !inside;
      }
    }
    return inside;
  };

  // all 4 must be inside
  return corners.every(inside);
}

/**
 * Extract only outer boundary loops (no holes) from a 2D triangulated mesh,
 * by throwing away any boundary loop that is contained in another.
 *
 * @param vertices Array of [x, y] coordinate pairs.
 * @param indices  Flat array of triangle indices; every consecutive 3 form one triangle.
 * @returns        Array of outer rings; each is an ordered [x, y][] tracing one CCW or CW boundary loop.
 */
export function extractOuterBoundaryRings(
  vertices: [number, number][],
  indices: ArrayLike<number>
): [number, number][][] {
  const makeKey = (i: number, j: number) => (i < j ? `${i}_${j}` : `${j}_${i}`);

  // 1) Count each undirected edge
  const edgeCount = new Map<string, number>();
  for (let t = 0; t < indices.length; t += 3) {
    const [i0, i1, i2] = [indices[t], indices[t + 1], indices[t + 2]];
    [
      [i0, i1],
      [i1, i2],
      [i2, i0],
    ].forEach(([a, b]) => {
      const key = makeKey(a, b);
      edgeCount.set(key, (edgeCount.get(key) || 0) + 1);
    });
  }

  // 2) Build adjacency for boundary edges (count === 1)
  const adj = new Map<number, number[]>();
  for (const [key, cnt] of edgeCount.entries()) {
    if (cnt !== 1) continue;
    const [a, b] = key.split("_").map(Number);
    if (!adj.has(a)) adj.set(a, []);
    if (!adj.has(b)) adj.set(b, []);
    adj.get(a)!.push(b);
    adj.get(b)!.push(a);
  }

  // 3) Traverse all boundary loops
  const visitedEdges = new Set<string>();
  const rings: [number, number][][] = [];

  for (const [start, neighbors] of adj.entries()) {
    for (const neighbor of neighbors) {
      const startKey = makeKey(start, neighbor);
      if (visitedEdges.has(startKey)) continue;

      const loopIdx: number[] = [start];
      let prev = start;
      let curr = neighbor;
      visitedEdges.add(startKey);
      loopIdx.push(curr);

      while (curr !== start) {
        const nbrs = adj.get(curr)!;
        let next: number | undefined;
        for (const w of nbrs) {
          const key = makeKey(curr, w);
          if (w !== prev && !visitedEdges.has(key)) {
            next = w;
            break;
          }
        }
        if (next === undefined) break; // dead end
        visitedEdges.add(makeKey(curr, next));
        loopIdx.push(next);
        prev = curr;
        curr = next;
      }

      rings.push(loopIdx.map((vi) => vertices[vi]));
    }
  }

  // 4) Point-in-polygon (ray-casting) test
  function pointInPoly(
    pt: [number, number],
    poly: [number, number][]
  ): boolean {
    const [px, py] = pt;
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const [xi, yi] = poly[i];
      const [xj, yj] = poly[j];
      // check if edge crosses the horizontal line at py
      const intersects =
        yi > py !== yj > py && px < ((xj - xi) * (py - yi)) / (yj - yi) + xi;
      if (intersects) inside = !inside;
    }
    return inside;
  }

  // 5) Filter out any ring that lies inside another ring
  return rings.filter((ring, i) => {
    // pick any point on this ring
    const testPt = ring[0];
    // if it's inside *any* other ring, it's a hole
    for (let j = 0; j < rings.length; j++) {
      if (i === j) continue;
      if (pointInPoly(testPt, rings[j])) {
        return false; // discard
      }
    }
    return true;
  });
}

/**
 * Check if a triangle mesh completely covers a given axis-aligned 2D bounding box,
 * reading each vertex as [x, y, z].
 *
 * @param vertices Flat number array [x0, y0, z0, x1, y1, z1, …, xN-1, yN-1, zN-1].
 * @param indices  Flat number array of triangle index triples [i0, i1, i2, …].
 * @param bbox     [xmin, ymin, xmax, ymax].
 * @param tol      Tolerance for floating-point comparisons (default 1e-9).
 * @returns        True if the mesh fully covers the bbox, false otherwise.
 */
function meshCoversBBox3D(
  vertices: ArrayLike<number>,
  indices: ArrayLike<number>,
  bbox: [number, number, number, number],
  tol = 1e-9
): boolean {
  const [xmin, ymin, xmax, ymax] = bbox;
  if (xmin >= xmax || ymin >= ymax) return false;
  if (indices.length < 3) return false;
  if (vertices.length % 3 !== 0)
    throw new Error("Vertices array must be a multiple of 3");

  // Collect all y-events: bbox edges + each vertex's y (clamped)
  const ySet = new Set<number>();
  ySet.add(ymin);
  ySet.add(ymax);
  for (let vi = 1; vi < vertices.length; vi += 3) {
    const y = vertices[vi];
    if (y > ymin + tol && y < ymax - tol) {
      ySet.add(y);
    }
  }
  const yList = Array.from(ySet).sort((a, b) => a - b);

  type Interval = { x0: number; x1: number };

  /** Test coverage on a horizontal line at yMid */
  function coversAtY(yMid: number): boolean {
    const intervals: Interval[] = [];

    for (let t = 0; t + 2 < indices.length; t += 3) {
      // triangle vertex indices
      const i0 = indices[t] * 3;
      const i1 = indices[t + 1] * 3;
      const i2 = indices[t + 2] * 3;

      const x0 = vertices[i0],
        y0 = vertices[i0 + 1];
      const x1 = vertices[i1],
        y1 = vertices[i1 + 1];
      const x2 = vertices[i2],
        y2 = vertices[i2 + 1];

      // skip if the scanline is outside the triangle's Y-range
      const minY = Math.min(y0, y1, y2);
      const maxY = Math.max(y0, y1, y2);
      if (yMid <= minY + tol || yMid >= maxY - tol) continue;

      // find intersections with triangle edges
      const xs: number[] = [];
      const pts: [number, number][] = [
        [x0, y0],
        [x1, y1],
        [x2, y2],
      ];
      for (let e = 0; e < 3; ++e) {
        const [xa, ya] = pts[e];
        const [xb, yb] = pts[(e + 1) % 3];
        if ((ya > yMid && yb < yMid) || (yb > yMid && ya < yMid)) {
          const t = (yMid - ya) / (yb - ya);
          xs.push(xa + t * (xb - xa));
        }
      }
      if (xs.length < 2) continue;
      xs.sort((a, b) => a - b);
      // clamp interval to bbox
      const left = Math.max(xs[0], xmin);
      const right = Math.min(xs[xs.length - 1], xmax);
      if (right > left + tol) {
        intervals.push({ x0: left, x1: right });
      }
    }

    if (intervals.length === 0) return false;
    // merge intervals
    intervals.sort((a, b) => a.x0 - b.x0);
    let curr = intervals[0];
    const merged: Interval[] = [curr];
    for (let i = 1; i < intervals.length; ++i) {
      const iv = intervals[i];
      if (iv.x0 <= curr.x1 + tol) {
        curr.x1 = Math.max(curr.x1, iv.x1);
      } else {
        curr = { x0: iv.x0, x1: iv.x1 };
        merged.push(curr);
      }
    }
    // check full coverage
    if (merged[0].x0 > xmin + tol) return false;
    if (merged[merged.length - 1].x1 < xmax - tol) return false;
    return true;
  }

  // test each slab between consecutive ys
  for (let i = 0; i + 1 < yList.length; ++i) {
    const yLo = yList[i],
      yHi = yList[i + 1];
    if (yHi - yLo < tol) continue;
    const yMid = 0.5 * (yLo + yHi);
    if (!coversAtY(yMid)) return false;
  }

  return true;
}

/**
 * Compute skirt (border) vertex indices for an axis-aligned rectangular mesh,
 * and split each border into connected segments following mesh triangles.
 *
 * @param vertices  Flat vertex array: [x, y, z, x, y, z, ...].
 * @param indices   Triangle index array: [i0, i1, i2, i0, i1, i2, ...].
 * @param stride    Number of components per vertex (default = 3).
 * @param eps       Tolerance for bounding comparisons (default = 1e-6).
 * @returns         Object containing raw border lists and connected segments.
 */
export interface SkirtIndices {
  east: number[];
  north: number[];
  west: number[];
  south: number[];
  eastSegments: number[][];
  northSegments: number[][];
  westSegments: number[][];
  southSegments: number[][];
}

export function getSkirtIndicesFlat(
  vertices: number[],
  indices: number[],
  stride = 3,
  eps = 1e-6
): SkirtIndices {
  const n = Math.floor(vertices.length / stride);
  if (n < 1) {
    return {
      east: [],
      north: [],
      west: [],
      south: [],
      eastSegments: [],
      northSegments: [],
      westSegments: [],
      southSegments: [],
    };
  }

  // 1) Find bounding extents in X and Y
  let minX = Infinity,
    maxX = -Infinity;
  let minY = Infinity,
    maxY = -Infinity;
  for (let i = 0; i < n; i++) {
    const x = vertices[i * stride];
    const y = vertices[i * stride + 1];
    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
    minY = Math.min(minY, y);
    maxY = Math.max(maxY, y);
  }

  // 2) Collect raw border vertex indices
  const east: number[] = [];
  const north: number[] = [];
  const west: number[] = [];
  const south: number[] = [];

  for (let i = 0; i < n; i++) {
    const x = vertices[i * stride];
    const y = vertices[i * stride + 1];
    if (Math.abs(x - maxX) <= eps) east.push(i);
    if (Math.abs(x - minX) <= eps) west.push(i);
    if (Math.abs(y - maxY) <= eps) north.push(i);
    if (Math.abs(y - minY) <= eps) south.push(i);
  }

  // 3) Sort each side for consistent ordering
  east.sort((a, b) => vertices[a * stride + 1] - vertices[b * stride + 1]);
  west.sort((a, b) => vertices[a * stride + 1] - vertices[b * stride + 1]);
  north.sort((a, b) => vertices[a * stride] - vertices[b * stride]);
  south.sort((a, b) => vertices[a * stride] - vertices[b * stride]);

  // 4) Build vertex adjacency from triangles
  const adj = new Map<number, Set<number>>();
  function link(u: number, v: number) {
    if (!adj.has(u)) adj.set(u, new Set());
    if (!adj.has(v)) adj.set(v, new Set());
    adj.get(u)!.add(v);
    adj.get(v)!.add(u);
  }
  for (let k = 0; k < indices.length; k += 3) {
    const [i0, i1, i2] = [indices[k], indices[k + 1], indices[k + 2]];
    link(i0, i1);
    link(i1, i2);
    link(i2, i0);
  }

  // 5) Split a sorted border into connected runs
  function splitConnected(sorted: number[]): number[][] {
    const segments: number[][] = [];
    if (sorted.length === 0) return segments;

    let current = [sorted[0]];
    for (let i = 1; i < sorted.length; i++) {
      const prev = sorted[i - 1];
      const cur = sorted[i];
      // if they share an edge, remain in same segment
      if (adj.get(prev)?.has(cur)) {
        current.push(cur);
      } else {
        segments.push(current);
        current = [cur];
      }
    }
    segments.push(current);
    return segments;
  }

  // 6) Compute connected segments for each side
  const eastSegments = splitConnected(east);
  const westSegments = splitConnected(west);
  const northSegments = splitConnected(north);
  const southSegments = splitConnected(south);

  return {
    east,
    north,
    west,
    south,
    eastSegments,
    northSegments,
    westSegments,
    southSegments,
  };
}

export function splitSharpEdges(creaseAngle: number = Math.PI / 3): Transform {
  return (document) => {
    for (const primitive of document
      .getRoot()
      .listMeshes()
      .flatMap((mesh) => mesh.listPrimitives())) {
      const indicesAccessor = primitive.getIndices();
      const positionAccessor = primitive.getAttribute("POSITION");

      if (!indicesAccessor || !positionAccessor) continue;

      const indices = indicesAccessor.getArray();
      const positions = positionAccessor.getArray();

      if (!indices || !positions) continue;

      const { indices: newIndices, vertices: newVertices } =
        splitVerticesAtSharpEdges(
          Array.from(positions),
          Array.from(indices),
          creaseAngle
        );

      indicesAccessor.setArray(new Uint32Array(newIndices));
      positionAccessor.setArray(new Float32Array(newVertices));

      compactPrimitive(primitive);
    }
  };
}

interface Mesh {
  vertices: number[]; // flat [x, y, z,  x, y, z, …]
  indices: number[]; // flat [i0, i1, i2,  i0, i1, i2, …]
}

/**
 * Splits vertices at sharp edges, computing face normals internally.
 *
 * @param vertices    Flat array of xyz positions.
 * @param indices     Flat array of triangle indices.
 * @param creaseAngle Maximum smoothing angle (in radians). Edges whose
 *                    adjacent face normals differ by more than this
 *                    angle will be “sharp” and cause a vertex split.
 */
function splitVerticesAtSharpEdges(
  vertices: number[],
  indices: number[],
  creaseAngle: number
): Mesh {
  const triCount = indices.length / 3;
  const faceNormals: number[] = new Array(triCount * 3);

  // 1) compute per-triangle normals
  for (let t = 0; t < triCount; t++) {
    const i0 = indices[3 * t + 0];
    const i1 = indices[3 * t + 1];
    const i2 = indices[3 * t + 2];
    const p0x = vertices[3 * i0 + 0],
      p0y = vertices[3 * i0 + 1],
      p0z = vertices[3 * i0 + 2];
    const p1x = vertices[3 * i1 + 0],
      p1y = vertices[3 * i1 + 1],
      p1z = vertices[3 * i1 + 2];
    const p2x = vertices[3 * i2 + 0],
      p2y = vertices[3 * i2 + 1],
      p2z = vertices[3 * i2 + 2];

    const e1x = p1x - p0x,
      e1y = p1y - p0y,
      e1z = p1z - p0z;
    const e2x = p2x - p0x,
      e2y = p2y - p0y,
      e2z = p2z - p0z;

    // cross(e1, e2)
    let nx = e1y * e2z - e1z * e2y;
    let ny = e1z * e2x - e1x * e2z;
    let nz = e1x * e2y - e1y * e2x;

    // normalize
    const len = Math.hypot(nx, ny, nz) || 1;
    nx /= len;
    ny /= len;
    nz /= len;

    faceNormals[3 * t + 0] = nx;
    faceNormals[3 * t + 1] = ny;
    faceNormals[3 * t + 2] = nz;
  }

  const cosThresh = Math.cos(creaseAngle);

  // 2) build edge → triangles map
  const edges = new Map<string, number[]>();
  for (let t = 0; t < triCount; t++) {
    const [i0, i1, i2] = [
      indices[3 * t + 0],
      indices[3 * t + 1],
      indices[3 * t + 2],
    ];
    for (const [a, b] of [
      [i0, i1],
      [i1, i2],
      [i2, i0],
    ] as [number, number][]) {
      const key = a < b ? `${a}_${b}` : `${b}_${a}`;
      const arr = edges.get(key);
      if (arr) arr.push(t);
      else edges.set(key, [t]);
    }
  }

  // 3) detect sharp edges
  const sharpEdges = new Set<string>();
  for (const [key, tris] of edges) {
    if (tris.length === 2) {
      const [tA, tB] = tris;
      const na = faceNormals.slice(3 * tA, 3 * tA + 3) as [
        number,
        number,
        number
      ];
      const nb = faceNormals.slice(3 * tB, 3 * tB + 3) as [
        number,
        number,
        number
      ];
      const dot = na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2];
      if (dot < cosThresh) sharpEdges.add(key);
    }
  }

  // 4) collect incident triangles per vertex
  const vertToTris = new Map<number, number[]>();
  for (let t = 0; t < triCount; t++) {
    for (let j = 0; j < 3; j++) {
      const v = indices[3 * t + j];
      const arr = vertToTris.get(v);
      if (arr) arr.push(t);
      else vertToTris.set(v, [t]);
    }
  }

  // 5) split vertices: build new vertex list and index mapping
  const cornerMap = new Map<string, number>();
  const newVerts: number[] = [];
  const newIndices = new Array(indices.length);

  for (const [v, tris] of vertToTris) {
    // adjacency of tris around this vertex via non-sharp edges
    const adj = new Map<number, Set<number>>();
    tris.forEach((t) => adj.set(t, new Set()));

    for (let i = 0; i < tris.length; i++) {
      for (let j = i + 1; j < tris.length; j++) {
        const tA = tris[i],
          tB = tris[j];
        const idxA = [
          indices[3 * tA],
          indices[3 * tA + 1],
          indices[3 * tA + 2],
        ];
        const idxB = [
          indices[3 * tB],
          indices[3 * tB + 1],
          indices[3 * tB + 2],
        ];
        if (!idxA.includes(v) || !idxB.includes(v)) continue;

        const othersA = idxA.filter((x) => x !== v);
        const common = othersA.find((x) => idxB.includes(x));
        if (common !== undefined) {
          const eKey = v < common ? `${v}_${common}` : `${common}_${v}`;
          if (!sharpEdges.has(eKey)) {
            adj.get(tA)!.add(tB);
            adj.get(tB)!.add(tA);
          }
        }
      }
    }

    // flood-fill smoothing groups
    const seen = new Set<number>();
    for (const start of tris) {
      if (seen.has(start)) continue;
      const stack = [start];
      seen.add(start);

      // gather this component
      const comp: number[] = [];
      while (stack.length) {
        const cur = stack.pop()!;
        comp.push(cur);
        for (const nb of adj.get(cur)!) {
          if (!seen.has(nb)) {
            seen.add(nb);
            stack.push(nb);
          }
        }
      }

      // create one new vertex for this smoothing group
      const newIndex = newVerts.length / 3;
      const px = vertices[3 * v],
        py = vertices[3 * v + 1],
        pz = vertices[3 * v + 2];
      newVerts.push(px, py, pz);

      // map each corner in the component to this new index
      for (const f of comp) {
        for (let corner = 0; corner < 3; corner++) {
          if (indices[3 * f + corner] === v) {
            cornerMap.set(`${f}_${corner}`, newIndex);
          }
        }
      }
    }
  }

  // 6) rebuild indices
  for (let t = 0; t < triCount; t++) {
    for (let j = 0; j < 3; j++) {
      const key = `${t}_${j}`;
      const ni = cornerMap.get(key);
      if (ni === undefined) {
        throw new Error(`Missing mapping for face ${t} corner ${j}`);
      }
      newIndices[3 * t + j] = ni;
    }
  }

  return { vertices: newVerts, indices: newIndices };
}

export function betterWeld(opts: { eps: number } = { eps: 0.1 }): Transform {
  return (document) => {
    for (const primitive of document
      .getRoot()
      .listMeshes()
      .flatMap((mesh) => mesh.listPrimitives())) {
      const positionAttribute = primitive.getAttribute("POSITION");
      const indexAttribute = primitive.getIndices();

      if (!positionAttribute || !indexAttribute) continue;

      const positions = positionAttribute.getArray();
      const indices = indexAttribute.getArray();

      if (!positions || !indices) continue;

      const { vertices: newPositions, indices: newIndices } = weld(
        positions,
        indices,
        opts.eps
      );

      positionAttribute.setArray(newPositions);
      indexAttribute.setArray(newIndices);
    }
  };
}

/**
 * Welds a mesh by merging vertices within an epsilon distance
 * and reorienting triangles so shared edges have opposite winding.
 */
export interface WeldResult {
  vertices: Float32Array;
  indices: Uint32Array;
}

/**
 * Weld operation on a mesh.
 * @param vertices - Flat array of 3D vertex positions [x0,y0,z0, x1,y1,z1, ...]
 * @param indices - Flat array of triangle indices [i0,j0,k0, i1,j1,k1, ...]
 * @param eps - Distance threshold for merging vertices
 * @returns Welded mesh with unique vertices and updated indices
 */
function weld(
  vertices: ArrayLike<number>,
  indices: ArrayLike<number>,
  eps: number = 1e-4
): WeldResult {
  const vertexCount = vertices.length / 3;
  const grid = new Map<string, number>();
  const unique: number[] = [];
  const remap = new Uint32Array(vertexCount);
  const eps2 = eps * eps;

  // Quantize position to a grid cell
  function key(x: number, y: number, z: number): string {
    const qx = Math.round(x / eps);
    const qy = Math.round(y / eps);
    const qz = Math.round(z / eps);
    return `${qx}_${qy}_${qz}`;
  }

  // Merge vertices
  for (let i = 0; i < vertexCount; i++) {
    const x = vertices[i * 3];
    const y = vertices[i * 3 + 1];
    const z = vertices[i * 3 + 2];
    const k = key(x, y, z);

    if (grid.has(k)) {
      const rep = grid.get(k)!;
      const rx = unique[rep * 3];
      const ry = unique[rep * 3 + 1];
      const rz = unique[rep * 3 + 2];
      const dx = x - rx;
      const dy = y - ry;
      const dz = z - rz;
      if (dx * dx + dy * dy + dz * dz <= eps2) {
        remap[i] = rep;
        continue;
      }
      // collision: will create new
    }

    const id = unique.length / 3;
    grid.set(k, id);
    unique.push(x, y, z);
    remap[i] = id;
  }

  // Build new index list, skipping degenerate tris
  const out: number[] = [];
  for (let f = 0; f < indices.length; f += 3) {
    const i0 = remap[indices[f] as number];
    const i1 = remap[indices[f + 1] as number];
    const i2 = remap[indices[f + 2] as number];
    if (i0 === i1 || i1 === i2 || i2 === i0) continue;
    out.push(i0, i1, i2);
  }

  // Reorient triangles for consistent winding along shared edges
  reorient(out);

  return {
    vertices: new Float32Array(unique),
    indices: new Uint32Array(out),
  };
}

/**
 * Ensures that for each pair of triangles sharing an edge,
 * their winding order along that edge is opposite.
 */
function reorient(indices: number[]): void {
  interface EdgeRec {
    face: number;
    dir: number;
  }
  const edgeMap = new Map<string, EdgeRec[]>();
  const faceCount = indices.length / 3;

  // Collect edges
  for (let f = 0; f < faceCount; f++) {
    const i = f * 3;
    const [a, b, c] = [indices[i], indices[i + 1], indices[i + 2]];
    const edges: [number, number][] = [
      [a, b],
      [b, c],
      [c, a],
    ];
    for (const [u, v] of edges) {
      const key = u < v ? `${u}_${v}` : `${v}_${u}`;
      const dir = u < v ? 1 : -1;
      if (!edgeMap.has(key)) edgeMap.set(key, []);
      edgeMap.get(key)!.push({ face: f, dir });
    }
  }

  // Flip triangles sharing edges with same direction
  for (const records of edgeMap.values()) {
    if (records.length === 2 && records[0].dir === records[1].dir) {
      const f = records[1].face;
      const idx = f * 3;
      // Swap 2nd and 3rd indices to flip winding
      [indices[idx + 1], indices[idx + 2]] = [
        indices[idx + 2],
        indices[idx + 1],
      ];
    }
  }
}
