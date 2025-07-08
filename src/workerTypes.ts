export type WorkerPayloads = WorkerInitPayload | WorkerWorkPayload;

export type WorkerInitPayload = {
  type: "init";
  tree: Record<string, any>;
  outputFolder: string;
  skipLevel: number;
};

export type WorkerWorkPayload = {
  type: "work";
  x: number;
  y: number;
  zoom: number;
};

export type WorkerWorkReturnType = {
  path: string;
  file: Buffer;
};
