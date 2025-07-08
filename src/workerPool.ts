import { Worker } from "worker_threads";

export class WorkerPool<
  Input extends Record<string, any>,
  Output
> {
  private readonly _workers: Worker[] = [];
  private readonly _busyList: boolean[] = [];
  private readonly _workQueue: {
    resolve: (value: any) => void;
    reject: (reason: any) => void;
    args: Record<string, any>;
  }[] = [];
  private _terminated: boolean;

  constructor(createWorker: () => Worker, workerCount: number = 16) {
    this._terminated = false;
    for (let i = 0; i < workerCount; i++) {
      this._workers.push(createWorker());
      this._busyList.push(false);
    }
  }

  async messageWorkers(args: Input): Promise<boolean> {
    const promises: Promise<any>[] = [];

    for (const worker of this._workers) {
      promises.push(
        new Promise((resolve, reject) => {
          worker.removeAllListeners();

          worker.on("message", (data) => {
            resolve(data);
          });
          worker.on("error", (error) => reject(error));
          worker.postMessage(args);
        })
      );
    }

    await Promise.all(promises);

    return true;
  }

  async run<A extends Input, B extends Output>(args: A): Promise<B> {
    if (this._terminated)
      throw new Error("WorkerPool was terminated, please use new WorkerPool!");

    let tempResolve: (value: any) => void = () => {};
    let tempReject: (reason: any) => void = () => {};
    const promise = new Promise<any>((resolve, reject) => {
      tempResolve = resolve;
      tempReject = reject;
    });
    this._workQueue.push({
      resolve: tempResolve,
      reject: tempReject,
      args,
    });
    this._scheduleWork();
    return promise;
  }

  private async _runWorker(
    index: number,
    args: Record<string, any>,
    resolve: (value: any) => void,
    reject: (value: any) => void
  ): Promise<any> {
    this._busyList[index] = true;

    const worker = this._workers[index];

    worker.removeAllListeners();

    worker.on("message", (data) => {
      resolve(data);
      this._busyList[index] = false;
      this._scheduleWork();
    });
    worker.on("error", (error) => {
      reject(error);
      this._busyList[index] = false;
      this._scheduleWork();
    });
    worker.postMessage(args);
  }

  private async _scheduleWork() {
    if (this._workQueue.length !== 0) {
      for (let i = 0; i < this._busyList.length; i++) {
        if (this._busyList[i]) continue;
        const { resolve, reject, args } = this._workQueue.pop()!;
        return this._runWorker(i, args, resolve, reject);
      }
    }
  }

  terminate() {
    for (const worker of this._workers) {
      worker.terminate();
    }
    this._terminated = true;
  }
}
