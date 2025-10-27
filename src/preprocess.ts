import { writeFile } from "fs/promises";
import gdal from "gdal-async";
import { basename, join } from "path";
import proj4 from "proj4";
import RBush from "rbush";
import { Box3, Vector3 } from "three";
import { triangulate3DPolygon } from "./helpers.js";
import { Tile, Tiles } from "./types.js";
import { processZip, yauzlOpen } from "./zip.js";
import { serializeAndWriteIndicesAndVertices } from "./customFormat.js";

export default async function preprocess(
  inputZip: string,
  outputFolder: string,
  progressCallback: (progress: number) => void,
  srcProj4?: string
) {
  type LineString = {
    type: "LineString";
    coordinates: [number, number, number][];
  };

  type Polygon = {
    type: "Polygon";
    coordinates: [number, number, number][][];
  };

  function arePointsEqual(
    a: [number, number, number],
    b: [number, number, number]
  ) {
    return a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
  }

  function buildGeometry(
    obj: LineString | Polygon | { type: string; coordinates: never }
  ) {
    switch (obj.type) {
      case "LineString":
        const firstPoint = obj.coordinates[0];
        const lastPoint = obj.coordinates.slice(-1)[0];

        if (arePointsEqual(firstPoint, lastPoint)) obj.coordinates.pop();

        return triangulate3DPolygon([obj.coordinates]);
      case "Polygon":
        for (const coordinates of obj.coordinates) {
          const firstPoint = coordinates[0];
          const lastPoint = coordinates.slice(-1)[0];

          if (arePointsEqual(firstPoint, lastPoint)) coordinates.pop();
        }

        return triangulate3DPolygon(obj.coordinates);
      default:
        console.error(`Unsupported GeoJSOn type found ${obj.type}`);

        return [];
    }
  }

  // Use Terrain db with triangles
  // Query over triangle set?



  const globalBoundingBox = new Box3();

  let index = 0;

  const tree = new RBush<Tile>();

  const zipFile = await yauzlOpen(inputZip, { lazyEntries: true });

  let total = zipFile.entryCount;

  await processZip(
    zipFile,
    // GDAL async supported 3d vector file formats
    [
      ".dxf",
      ".dgn",
      ".shp",
      ".geojson",
      ".gpkg",
      ".gpx",
      ".kml",
      ".kmz",
      ".mvt",
      ".pbf",
      ".sqlite",
      ".svg",
    ],
    async (fileName, data) => {
      try {
        const outputPath = join(outputFolder, basename(fileName) + ".bin");
        const set = gdal.open(data);

        const srcSRS = set.srs?.toProj4() ?? srcProj4;
        const destSRS3857 =
          "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs";
        if (!srcSRS)
          throw new Error("No srs found, please provide a srcProj4 string.");
        const transform3857 = proj4(srcSRS, destSRS3857);

        let walkingIndex = 0;

        let meshIndices: number[] = [];
        let meshVertices: number[] = [];

        let gt: Vector3 | null = null;
        const localBoundingBox = new Box3();

        for (const layer of set.layers) {
          layer.features.forEach((feature) => {
            const obj = feature.getGeometry().toObject() as
              | LineString
              | Polygon;

            const indices = buildGeometry(obj);

            const transformedCoordinates3857 = (
              obj.type === "Polygon"
                ? obj.coordinates.flatMap((d) => d)
                : obj.coordinates
            ).map((coordinates) => {
              return transform3857.forward(coordinates);
            });

            if (!indices) return;

            if (!gt) {
              gt = new Vector3(
                transformedCoordinates3857[0][0],
                transformedCoordinates3857[0][1],
                transformedCoordinates3857[0][2]
              );
            }

            const tempBox = new Box3();

            tempBox.setFromArray(transformedCoordinates3857.flat());

            localBoundingBox.union(tempBox);

            meshIndices.push(...indices.map((index) => index + walkingIndex));
            meshVertices.push(
              ...transformedCoordinates3857.flatMap((number) => [
                number[0] - gt!.x,
                number[1] - gt!.y,
                number[2] - gt!.z,
              ])
            );

            walkingIndex += transformedCoordinates3857.length;
          });
        }

        globalBoundingBox.union(localBoundingBox);

        tree.insert({
          minX: localBoundingBox.min.x,
          minY: localBoundingBox.min.y,
          minZ: localBoundingBox.min.z,
          maxX: localBoundingBox.max.x,
          maxY: localBoundingBox.max.y,
          maxZ: localBoundingBox.max.z,
          path: outputPath,
        });

        await serializeAndWriteIndicesAndVertices(
          outputPath,
          gt!,
          meshIndices,
          meshVertices
        );

        try {
          await writeFile(
            join(outputFolder, "metadata.json"),
            JSON.stringify({
              minX: globalBoundingBox.min.x,
              minY: globalBoundingBox.min.y,
              minZ: globalBoundingBox.min.z,
              maxX: globalBoundingBox.max.x,
              maxY: globalBoundingBox.max.y,
              maxZ: globalBoundingBox.max.z,
              tiles: tree.toJSON(),
            } as Tiles)
          );
        } catch (e) {
          console.error(e);
        }

        index++;

        progressCallback(index / total);
      } catch (e) {
        console.error(e);
      }
    }
  );
}
