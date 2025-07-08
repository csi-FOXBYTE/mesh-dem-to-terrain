import RBush from "rbush";

export type Tile = {
    path: string;
    minX: number;
    minY: number;
    minZ: number;
    maxX: number;
    maxY: number;
    maxZ: number;
  };
  
  export type Tiles = {
    minX: number;
    minY: number;
    minZ: number;
    maxX: number;
    maxY: number;
    maxZ: number;
    tiles: ReturnType<RBush<Tile>["toJSON"]>;
  };
  