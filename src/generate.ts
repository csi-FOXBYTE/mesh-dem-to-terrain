import { PromisePool } from "@supercharge/promise-pool";
import fs, { readFile } from "fs/promises";
import path, { dirname } from "path";
import proj4 from "proj4";
import slippyGrid from "slippy-grid";
import { Worker } from "worker_threads";
import { Tiles } from "./types.js";
import { WorkerPool } from "./workerPool.js";
import {
  WorkerInitPayload,
  WorkerWorkPayload,
  WorkerWorkReturnType,
} from "./workerTypes.js";

export default async function generate(
  outputFolder: string,
  progressCallback: (progress: number) => void,
  opts: {
    writeFile?: (
      path: string,
      file: Buffer | string,
      terrainTileInfos?: { x: number; y: number; zoom: number }
    ) => Promise<void>;
    threadCount?: number;
    endZoom?: number;
    startZoom?: number;
    skipLevel?: number;
  }
) {
  const {
    endZoom = 16,
    startZoom = 0,
    writeFile = async (path, file) => {
      await fs.mkdir(dirname(path), { recursive: true });
      await fs.writeFile(path, file);
    },
    threadCount = 4,
    skipLevel = 10,
  } = opts;

  const metadata = JSON.parse(
    (await readFile(path.join(path.resolve(outputFolder), "metadata.json"))).toString()
  ) as Tiles;

  const srcSRS4326 = "+proj=longlat +datum=WGS84 +no_defs +type=crs";
  const destSRS3857 =
    "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs";
  const transform3857To4326 = proj4(destSRS3857, srcSRS4326);

  const [gx0, gy0] = transform3857To4326.forward([
    metadata.minX,
    metadata.minY,
  ]);
  const [gx1, gy1] = transform3857To4326.forward([
    metadata.maxX,
    metadata.maxY,
  ]);

  const extents = [gx0, gy0, gx1, gy1] as const;
  // const extents = [9.9278,53.5364,9.9308,53.5378] as const;

  const grid = slippyGrid.all(extents, startZoom, endZoom);

  const levels = slippyGrid.levels(extents, startZoom, endZoom);

  const total = levels.reduce((prev, [xs, ys]) => {
    return prev + xs.length * ys.length;
  }, 0);

  await writeFile(
    path.join(outputFolder, "layer.json"),
    JSON.stringify(
      {
        tiles: ["{z}/{x}/{y}.terrain"],
        projection: "EPSG:3857",
        bounds: [-180, -90, 180, 90],
        format: "quantized-mesh-1.0",
        available: levels.map(([xs, ys]) => [
          {
            startX: xs[0],
            startY: ys[0],
            endX: xs.slice(-1)[0],
            endY: ys.slice(-1)[0],
          },
        ]),
      },
      undefined,
      4
    )
  );

  const workerPool = new WorkerPool(
    () => new Worker(new URL("./worker.js", import.meta.url), {}),
    threadCount
  );

  await workerPool.messageWorkers({
    type: "init",
    tree: metadata.tiles,
    outputFolder,
    skipLevel,
  } as WorkerInitPayload);

  let index = 0;

  await new PromisePool(grid)
    .withConcurrency(threadCount)
    .process(async ([x, y, zoom]) => {
      const { file, path } = await workerPool.run<
        WorkerWorkPayload,
        WorkerWorkReturnType
      >({
        type: "work",
        x,
        y,
        zoom,
      });
      await writeFile(path, file, { x, y, zoom });
      index++;
      progressCallback(index / total);
    });
}
