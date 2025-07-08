import type { Options } from "tsup";

export const tsup: Options = {
  splitting: true,
  sourcemap: true,
  platform: "node",
  minify: true,
  dts: {
    entry: "src/index.ts",
  },
  target: "es2020",
  format: ["esm"],
  bundle: false,
  clean: true,
  treeshake: true,
  entry: ["src/**/*.ts"],
  noExternal: [],
};
