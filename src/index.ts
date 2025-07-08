import generate from "./generate.js";
import preprocess from "./preprocess.js";

export { generate, preprocess };

(async () => {
  await preprocess("D:\\dxf.zip", "D:\\dxf_test", console.log, "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs");

  await generate("D:\\dxf_test", (progress) =>
    console.log((progress * 100).toFixed(2), "%")
  );

  process.exit(0);
})();
