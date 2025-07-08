# mesh-dem-to-terrain 🚀🥾

A Node.js library for converting DEM meshes (DXF, GeoJSON, etc.) into Cesium Terrain Tiles (quantized-mesh or heightmap), with flexible zoom-level control.

## Table of Contents

* [Features](#features)
* [Installation](#installation)
* [Usage](#usage)
* [API](#api)
* [Options Overview](#options-overview)
* [CLI Wrapper Example](#cli-wrapper-example)
* [Contributing](#contributing)
* [License](#license)

## 🎉 Features

* **🔄 Preprocess Mesh Archive**: Unpacks and prepares DXF/GeoJSON mesh archives into a working directory, reprojects coordinates using Proj4 strings. 🗺️
* **🗻 Terrain Tile Generation**: Generates Cesium-compatible terrain tiles (quantized-mesh or heightmap) across specified zoom levels. 📈
* **📐 Projection Support**: Accept any Proj4-formatted CRS string for flexible coordinate transformations. ✨
* **📊 Progress Callbacks**: Receive real-time progress updates (percentage or log messages). 📢
* **🧵 Multi-threading**: Customize worker threads for CPU-bound tasks. ⚙️

## 📥 Installation

```bash
npm install mesh-dem-to-terrain
```

## 💻 Usage

```js
import { preprocess, generate } from "mesh-dem-to-terrain/index.js";

(async () => {
  // Step 1: Unpack and reproject your mesh archive 🛠️
  await preprocess(
    "D:\\dxf.zip",                      // Path to DXF/GeoJSON archive
    "D:\\dxf_test",                    // Working directory for intermediate files
    console.log,                       // Progress callback (log messages)
    "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"  // Proj4 CRS string
  );

  // Step 2: Generate Cesium Terrain Tiles 🌄
  await generate(
    "D:\\dxf_test",                    // Working directory from preprocess
    (progress) => console.log((progress * 100).toFixed(2), "%"),  // Progress callback (percentage)
    { endZoom: 16, startZoom: 0, threadCount: 2 }  // Options: zoom range & threads
  );
})();
```

## ⚙️ API

### `preprocess(inputArchive, workDir, progressCallback, proj4String)`

* **inputArchive** `(string)` – Path to a ZIP archive containing DXF, GeoJSON, or other mesh files. 📂
* **workDir** `(string)` – Directory to extract and preprocess mesh data. 📁
* **progressCallback** `(function)` – Called with log messages or progress updates. 📢
* **proj4String** `(string)` – Proj4-formatted CRS string for coordinate transformation. 🌍

**Returns:** A promise that resolves when preprocessing is complete. ✅

### `generate(workDir, progressCallback, options)`

* **workDir** `(string)` – Directory prepared by `preprocess()`. 📂
* **progressCallback** `(function)` – Called with a number in \[0,1] indicating progress. 📈
* **options** `(object)`:

  * `startZoom` `(number)` – Minimum zoom level to generate (default: 0). ↗️
  * `endZoom` `(number)` – Maximum zoom level to generate (required). ↘️
  * `threadCount` `(number)` – Number of worker threads (default: number of CPU cores). 🧵

**Returns:** A promise that resolves when terrain tile generation is complete. ✅

## 🛠️ Options Overview

| Option        | Default            | Description                                |
| ------------- | ------------------ | ------------------------------------------ |
| `startZoom`   | `0`                | Minimum zoom level for generated tiles. ↗️ |
| `endZoom`     | *required*         | Maximum zoom level for generated tiles. ↘️ |
| `threadCount` | `os.cpus().length` | Number of parallel worker threads. 🧵      |

## 📜 CLI Wrapper Example

Wrap the library in a convenient CLI script:

```js
#!/usr/bin/env node
import path from "path";
import { preprocess, generate } from "mesh-dem-to-terrain/index.js";

const [,, archive, outDir, proj, startZ, endZ, threads] = process.argv;

(async () => {
  // Preprocess mesh archive 🛠️
  await preprocess(
    path.resolve(archive),
    path.resolve(outDir),
    console.log,
    proj
  );

  // Generate terrain tiles 🌄
  await generate(
    path.resolve(outDir),
    (p) => console.log(`${(p*100).toFixed(2)} %`),
    {
      startZoom: Number(startZ) || 0,
      endZoom:   Number(endZ) || 16,
      threadCount: Number(threads) || 4
    }
  );
})();
```

## 🤝 Contributing

Contributions are welcome! Please open issues or pull requests on the GitHub repository. 🙌

## 📄 License

MIT License. See [LICENSE](LICENSE) for details. 📝
