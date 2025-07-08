import { Cartesian3, Ellipsoid, WebMercatorProjection } from "cesium";
import _ from "lodash";
import { Box3, Sphere, Vector3 } from "three";
import { getNewellsNormal } from "./helpers.js";

const HEADER_SIZE = 3 * 8 + 2 * 4 + 4 * 8 + 3 * 8;

function zigZagEncode(value: number): number {
  // In JS, bitwise ops handle 32-bit signed integers.
  return (value << 1) ^ (value >> 31);
}

function writeArrayDeltaZigZag(
  buffer: Buffer,
  index: number,
  array: Uint16Array
) {
  let prev = 0;
  for (const element of array) {
    const current = element;
    const delta = current - prev; // delta from previous
    prev = current; // track for next iteration

    const zz = zigZagEncode(delta);
    buffer.writeUInt16LE(zz, index);
    index += 2;
  }

  return index;
}

function writeHighWaterMarkIndices(
  buffer: Buffer,
  index: number,
  indices: number[],
  encoded: boolean,
  needsBigIndices: boolean
) {
  const writeFunction = needsBigIndices
    ? (n: number, offset: number) => buffer.writeUInt32LE(n, offset)
    : (n: number, offset: number) => buffer.writeUInt16LE(n, offset);
  const increment = needsBigIndices ? 4 : 2;

  if (!encoded) {
    for (const i of indices) {
      writeFunction(i, index);
      index += increment;
    }

    return index;
  } else {
    let highest = 0;

    for (const i of indices) {
      const code = highest - i;

      writeFunction(code, index);

      index += increment;

      if (code === 0) {
        highest++;
      }
    }

    return index;
  }
}

export class QuantizedMesh {
  private centerX: number;
  private centerY: number;
  private centerZ: number;

  private minimumHeight: number;
  private maximumHeight: number;

  private boundingSphereCenterX: number;
  private boundingSphereCenterY: number;
  private boundingSphereCenterZ: number;
  private boundingSphereRadius: number;

  private horizonOcclusionPointX: number;
  private horizonOcclusionPointY: number;
  private horizonOcclusionPointZ: number;

  private vertexCount: number;
  private u: Uint16Array = new Uint16Array();
  private v: Uint16Array = new Uint16Array();
  private height: Uint16Array = new Uint16Array();

  private triangleCount: number;
  private indices: number[];
  private eastIndices: number[];
  private westIndices: number[];
  private southIndices: number[];
  private northIndices: number[];

  private octEncodedNormalsCount: number = 0;
  private octEncodedNormals: number[] = [];

  private _computeHorizonOcclusionPoint(points: Vector3[]) {
    const cartesian = new Cartesian3();

    function computeMagnitude(
      point: Cartesian3,
      scaledSpaceDirectionToPoint: Vector3
    ) {
      const p = Ellipsoid.WGS84.transformPositionToScaledSpace(point);
      const scaledSpacePoint = new Vector3(p.x, p.y, p.z);
      let magnitudeSquared = scaledSpacePoint.lengthSq();
      let magnitude = Math.sqrt(magnitudeSquared);
      const direction = scaledSpacePoint.divideScalar(magnitude);

      magnitudeSquared = Math.max(1.0, magnitudeSquared);
      magnitude = Math.max(1.0, magnitude);

      const cosAlpha = direction.dot(scaledSpaceDirectionToPoint);
      const sinAlpha = direction.cross(scaledSpaceDirectionToPoint).length();
      const cosBeta = 1.0 / magnitude;
      const sinBeta = Math.sqrt(magnitudeSquared - 1.0) * cosBeta;

      return 1.0 / (cosAlpha * cosBeta - sinAlpha * sinBeta);
    }

    const scaledSpaceDirectionToPointCartesian =
      Ellipsoid.WGS84.transformPositionToScaledSpace(
        new Cartesian3(this.centerX, this.centerY, this.centerZ)
      );
    const scaledSpaceDirectionToPoint = new Vector3(
      scaledSpaceDirectionToPointCartesian.x,
      scaledSpaceDirectionToPointCartesian.y,
      scaledSpaceDirectionToPointCartesian.z
    );

    const magnitudes: number[] = [];
    let maxMagnitude = -Infinity;

    for (const point of points) {
      cartesian.x = point.x;
      cartesian.y = point.y;
      cartesian.z = point.z;
      const magnitude = computeMagnitude(
        cartesian,
        scaledSpaceDirectionToPoint
      );

      maxMagnitude = Math.max(magnitude, maxMagnitude);
      magnitudes.push(magnitude);
    }

    scaledSpaceDirectionToPoint.multiplyScalar(maxMagnitude);

    const horizonOcclusionPoint =
      Ellipsoid.WGS84.transformPositionFromScaledSpace(
        new Cartesian3(
          scaledSpaceDirectionToPoint.x,
          scaledSpaceDirectionToPoint.y,
          scaledSpaceDirectionToPoint.z
        )
      );

    this.horizonOcclusionPointX = horizonOcclusionPoint.x;
    this.horizonOcclusionPointY = horizonOcclusionPoint.y;
    this.horizonOcclusionPointZ = horizonOcclusionPoint.z;
  }

  private _computeVertexNormals(points: Vector3[], indices: number[]) {
    const normals: Vector3[] = points.map(() => new Vector3(0, 0, 0));

    for (let i = 0; i < indices.length; ) {
      const index1 = indices[i++];
      const a = points[index1];
      const index2 = indices[i++];
      const b = points[index2];
      const index3 = indices[i++];
      const c = points[index3];

      const n = getNewellsNormal([a, b, c]);

      normals[index1].add(n);
      normals[index2].add(n);
      normals[index3].add(n);
    }

    normals.forEach((n) => {
      if (n.length() < 1e-6) {
        n.set(0, 0, 1);
      } else {
        n.normalize();
      }
    }); 

    const octEncodedNormals: number[] = [];

    const toSnorm = (v: number) =>
      Math.round((_.clamp(v, -1.0, 1.0) * 0.5 + 0.5) * 255.0);

    const sign = (v: number) => (v < 0 ? -1 : 1);

    for (const normal of normals) {
      const l1Norm =
        Math.abs(normal.x) + Math.abs(normal.y) + Math.abs(normal.z);

      let px = normal.x / l1Norm;
      let py = normal.y / l1Norm;

      if (normal.z < 0) {
        const rx = 1 - Math.abs(py) * sign(px);
        const ry = 1 - Math.abs(px) * sign(py);
        px = rx;
        py = ry;
      }

      octEncodedNormals.push(toSnorm(px), toSnorm(py));
    }

    this.octEncodedNormals = octEncodedNormals;
    this.octEncodedNormalsCount = octEncodedNormals.length;
  }

  private createPositionVectors(vertices: number[], bbox: Box3) {
    const us = new Uint16Array(vertices.length / 3);
    const vs = new Uint16Array(vertices.length / 3);
    const heights = new Uint16Array(vertices.length / 3);

    let walkingIndex = 0;

    for (let i = 0; i < vertices.length; ) {
      const x = vertices[i++];
      const y = vertices[i++];
      const heightUnnormalized = vertices[i++];

      const u = Math.round(
        ((x - bbox.min.x) / Math.max(bbox.max.x - bbox.min.x, 0.001)) * 32767
      );
      const v = Math.round(
        ((y - bbox.min.y) / Math.max(bbox.max.y - bbox.min.y, 0.001)) * 32767
      );
      const height = Math.round(
        ((heightUnnormalized - this.minimumHeight) /
          Math.max(this.maximumHeight - this.minimumHeight, 0.001)) *
          32767
      );

      us[walkingIndex] = u;
      vs[walkingIndex] = v;
      heights[walkingIndex] = height;

      walkingIndex++;
    }

    this.u = us;
    this.v = vs;
    this.height = heights;
  }

  constructor(
    vertices: number[],
    indices: number[],
    eastIndices: number[],
    westIndices: number[],
    southIndices: number[],
    northIndices: number[],
    bbox3857: Box3
  ) {
    this.minimumHeight = Infinity;
    this.maximumHeight = -Infinity;

    const bbox4978 = new Box3();

    const points: Vector3[] = []; // points in 4978 (ECEF)

    for (let i = 0; i < vertices.length; ) {
      const vert = [vertices[i++], vertices[i++], vertices[i++]] as [
        number,
        number,
        number
      ];

      const cartographic = new WebMercatorProjection(Ellipsoid.WGS84).unproject(
        new Cartesian3(...vert)
      );

      const cartesian = Ellipsoid.WGS84.cartographicToCartesian(cartographic);

      points.push(new Vector3(cartesian.x, cartesian.y, cartesian.z));

      bbox4978.expandByPoint(points[points.length - 1]);

      this.minimumHeight = Math.min(this.minimumHeight, vert[2]);
      this.maximumHeight = Math.max(this.maximumHeight, vert[2]);
    }

    bbox3857.max.z = this.maximumHeight;
    bbox3857.min.z = this.minimumHeight;

    const bbsphere = new Sphere().setFromPoints(points);
    const center = new Vector3();

    bbox4978.getCenter(center);

    this.indices = indices;
    this.eastIndices = eastIndices;
    this.westIndices = westIndices;
    this.southIndices = southIndices;
    this.northIndices = northIndices;

    this.centerX = center.x;
    this.centerY = center.y;
    this.centerZ = center.z;

    this.boundingSphereCenterX = bbsphere.center.x;
    this.boundingSphereCenterY = bbsphere.center.y;
    this.boundingSphereCenterZ = bbsphere.center.z;
    this.boundingSphereRadius = bbsphere.radius;

    this.horizonOcclusionPointX = NaN;
    this.horizonOcclusionPointY = -Infinity;
    this.horizonOcclusionPointZ = NaN;

    this.vertexCount = vertices.length / 3;
    this.triangleCount = indices.length / 3;

    this.createPositionVectors(vertices, bbox3857);

    this._computeVertexNormals(points, indices);

    this._computeHorizonOcclusionPoint(points);
  }

  async serialize() {
    const NEEDS_BIG_INDICES = this.triangleCount > 2 ** 16;

    const VERTEX_SIZE =
      4 + // vertex count
      2 * this.vertexCount + // us
      2 * this.vertexCount + // vs
      2 * this.vertexCount; // heights;

    const INDEX_SIZE = NEEDS_BIG_INDICES ? 4 : 2;

    const PADDING = (HEADER_SIZE + VERTEX_SIZE) % INDEX_SIZE;

    const totalBytes =
      HEADER_SIZE +
      VERTEX_SIZE +
      PADDING +
      4 + // triangle count
      INDEX_SIZE * this.indices.length + // indices
      4 + // west vertex count
      INDEX_SIZE * this.westIndices.length + // west indices
      4 + // south vertex count
      INDEX_SIZE * this.southIndices.length + // south indices
      4 + // east vertex count
      INDEX_SIZE * this.eastIndices.length + // east indices
      4 + // noth vertex count
      INDEX_SIZE * this.northIndices.length + // north indices
      1 + // normals extension id
      4 + // normals length
      1 * this.vertexCount * 2; // oct encoded normals

    const outBuffer = Buffer.alloc(totalBytes);
    let index = 0;

    // center components (3 x double)

    outBuffer.writeDoubleLE(this.centerX, index);
    index += 8;

    outBuffer.writeDoubleLE(this.centerY, index);
    index += 8;

    outBuffer.writeDoubleLE(this.centerZ, index);
    index += 8;

    // min, max height (2 x float)
    outBuffer.writeFloatLE(this.minimumHeight, index);
    index += 4;

    outBuffer.writeFloatLE(this.maximumHeight, index);
    index += 4;

    // bounding sphere center (3 x double)
    outBuffer.writeDoubleLE(this.boundingSphereCenterX, index);
    index += 8;

    outBuffer.writeDoubleLE(this.boundingSphereCenterY, index);
    index += 8;

    outBuffer.writeDoubleLE(this.boundingSphereCenterZ, index);
    index += 8;

    // bounding sphere radius (double)
    outBuffer.writeDoubleLE(this.boundingSphereRadius, index);
    index += 8;

    // horizon occlusion point (3 x double)
    outBuffer.writeDoubleLE(this.horizonOcclusionPointX, index);
    index += 8;

    outBuffer.writeDoubleLE(this.horizonOcclusionPointY, index);
    index += 8;

    outBuffer.writeDoubleLE(this.horizonOcclusionPointZ, index);
    index += 8;

    // vertex count (uint32)
    outBuffer.writeUInt32LE(this.vertexCount, index);
    index += 4;

    // 8. Vertex data (u, v, heights) each as int16, with correct offset
    index = writeArrayDeltaZigZag(outBuffer, index, this.u);
    index = writeArrayDeltaZigZag(outBuffer, index, this.v);
    index = writeArrayDeltaZigZag(outBuffer, index, this.height);

    // Add padding pre indices to ensure correct byte alignment
    for (let i = 0; i < PADDING; i++) {
      outBuffer.writeUint8(0, index++);
    }

    if (index % INDEX_SIZE !== 0) throw new Error("Needs padding!");

    // 9. Index count and index data
    outBuffer.writeUint32LE(this.triangleCount, index);
    index += 4;

    index = writeHighWaterMarkIndices(
      outBuffer,
      index,
      this.indices,
      true,
      NEEDS_BIG_INDICES
    );

    // 9. West Index count and index data (example)
    // Typically these could be 16-bit or 32-bit.
    // We'll assume they're small enough for 16-bit:
    outBuffer.writeUint32LE(this.westIndices.length, index);
    index += 4;

    index = writeHighWaterMarkIndices(
      outBuffer,
      index,
      this.westIndices,
      false,
      NEEDS_BIG_INDICES
    );

    // 9. South Index count and index data (example)
    // Typically these could be 16-bit or 32-bit.
    // We'll assume they're small enough for 16-bit:
    outBuffer.writeUint32LE(this.southIndices.length, index);
    index += 4;

    index = writeHighWaterMarkIndices(
      outBuffer,
      index,
      this.southIndices,
      false,
      NEEDS_BIG_INDICES
    );

    // 9. East Index count and index data (example)
    // Typically these could be 16-bit or 32-bit.
    // We'll assume they're small enough for 16-bit:
    outBuffer.writeUint32LE(this.eastIndices.length, index);
    index += 4;

    index = writeHighWaterMarkIndices(
      outBuffer,
      index,
      this.eastIndices,
      false,
      NEEDS_BIG_INDICES
    );

    // 9. North Index count and index data (example)
    // Typically these could be 16-bit or 32-bit.
    // We'll assume they're small enough for 16-bit:
    outBuffer.writeUint32LE(this.northIndices.length, index);
    index += 4;

    index = writeHighWaterMarkIndices(
      outBuffer,
      index,
      this.northIndices,
      false,
      NEEDS_BIG_INDICES
    );

    // Normals
    outBuffer.writeUint8(1, index++);
    outBuffer.writeUint32LE(this.octEncodedNormalsCount * 2, index);
    index += 4;
    for (const element of this.octEncodedNormals) {
      outBuffer.writeUint8(element, index++);
    }

    if (index !== totalBytes) throw new Error("Buffer length didnt match!");

    return outBuffer;
  }
}
