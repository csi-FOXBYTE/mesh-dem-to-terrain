import { readFile, writeFile } from "fs/promises";
import { Vector3 } from "three";

export async function readAndDeserializeIndicesAndVertices(input: string) {
    const file = await readFile(input);

    const m = file.readUint8(0);
    const e = file.readUint8(1);
    const s = file.readUint8(2);
    const h = file.readUint8(3);

    if (m !== 0x4d || e !==  0x45 || s !== 0x53 || h !== 0x48) {
        throw new Error("Read file is not of needed format!");
    }

    const gt = new Vector3();

    const x = file.readDoubleLE(4);
    const y = file.readDoubleLE(12);
    const z = file.readDoubleLE(20);

    const indexStart = 28;

    const indicesLength = file.readUInt32LE(indexStart);

    const indices = new Uint32Array(indicesLength);

    const isBigIndex = indicesLength > 2 ** 16;
    const indexStep = isBigIndex ? 4 : 2;

    if (isBigIndex) {
        for (let i = 0; i < indicesLength; i++) {
            indices[i] = file.readUint32LE(indexStart + 4 + i * indexStep);
        }
    } else {
        for (let i = 0; i < indicesLength; i++) {
            indices[i] = file.readUint16LE(indexStart + 4 + i * indexStep);
        }
    }

    const vertexStart = indexStart + 4 + indicesLength * indexStep;

    const verticesLength = file.readUInt32LE(vertexStart);

    const vertices = new Float32Array(verticesLength);

    for (let i = 0; i < verticesLength; i++) {
        vertices[i] = file.readDoubleLE(vertexStart + 4 + i * 8);
    }

    gt.set(x, y, z);

    return { gt, indices, vertices };
}

export async function serializeAndWriteIndicesAndVertices(
    out: string,
    gt: Vector3,
    indices: number[],
    vertices: number[]
) {
    const headerBuffer = Buffer.from({ length: 4 });

    headerBuffer.writeUInt8(0x4d, 0); // M
    headerBuffer.writeUInt8(0x45, 1); // E
    headerBuffer.writeUInt8(0x53, 2); // S
    headerBuffer.writeUInt8(0x48, 3); // H

    const globalTransformBuffer = Buffer.from({ length: 24 });

    globalTransformBuffer.writeDoubleLE(gt.x, 0);
    globalTransformBuffer.writeDoubleLE(gt.y, 8);
    globalTransformBuffer.writeDoubleLE(gt.z, 16);

    const isBigIndex = !(indices.length <= 2 ** 16);
    const indicesStride = isBigIndex ? 4 : 2;

    const indicesBuffer = Buffer.from({ length: 4 + indices.length * indicesStride });

    indicesBuffer.writeUInt32LE(indices.length, 0);

    let indicesOffset = 4;

    if (isBigIndex) {
        for (const index of indices) {
            indicesBuffer.writeUint32LE(index, indicesOffset);
            indicesOffset += indicesStride;
        }
    } else {
        for (const index of indices) {
            indicesBuffer.writeUint16LE(index, indicesOffset);
            indicesOffset += indicesStride;
        }
    };

    if (indicesBuffer.length !== indicesOffset) throw new Error("Range mismatch for indices buffer!");


    const verticesBuffer = Buffer.from({ length: 4 + vertices.length * 8 });

    verticesBuffer.writeUint32LE(vertices.length, 0);

    let verticesOffset = 4;

    for (const vertex of vertices) {
        verticesBuffer.writeDoubleLE(vertex, verticesOffset);
        verticesOffset += 8;
    }

    if (verticesBuffer.length !== verticesOffset) throw new Error("Range mismatch for vertices buffer!");

    await writeFile(out, [headerBuffer, globalTransformBuffer, indicesBuffer, verticesBuffer])
}