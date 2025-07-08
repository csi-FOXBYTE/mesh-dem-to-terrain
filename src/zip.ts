import yauzl from "yauzl";

export const yauzlOpen = (path: string, opts: yauzl.Options) => {
  return new Promise<yauzl.ZipFile>((resolve, reject) =>
    yauzl.open(path, opts, (err, zipFile) => {
      resolve(zipFile);
      reject(err);
    })
  );
};

export async function processZip(
  zipFile: yauzl.ZipFile,
  suffixes: string[],
  onEntryBuffer: (fileName: string, data: Buffer) => Promise<void>
): Promise<void> {
  // Open the zip in lazy mode so we can control when entries are read

  return new Promise((resolve, reject) => {
    zipFile.readEntry(); // start reading

    zipFile.on("entry", (entry: yauzl.Entry) => {
      // Skip directories
      if (/\/$/.test(entry.fileName)) {
        return zipFile.readEntry();
      }

      // If it matches our suffix, extract to buffer
      if (suffixes.some((suffix) => entry.fileName.endsWith(suffix))) {
        zipFile.openReadStream(entry, (err, readStream) => {
          if (err) return reject(err);
          const chunks: Buffer[] = [];
          readStream!.on("data", (chunk: Buffer) => chunks.push(chunk));
          readStream!.on("end", async () => {
            const fileBuffer = Buffer.concat(chunks);
            try {
              zipFile.readEntry();
              await onEntryBuffer(entry.fileName, fileBuffer);
            } catch (cbErr) {
              reject(cbErr);
            }
          });
          readStream!.on("error", reject);
        });
      } else {
        // Not our file, skip it
        zipFile.readEntry();
      }
    });

    zipFile.on("end", () => {
      zipFile.close();
      resolve();
    });
    zipFile.on("error", reject);
  });
}
