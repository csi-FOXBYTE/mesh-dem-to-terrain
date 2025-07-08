import { Box3, Matrix4, Quaternion, Vector3, Vector3Like } from "three";
import * as earcut from "earcut";
import _ from "lodash";
import proj4 from "proj4";
import proj4List from "proj4-list";

const cesiumCartographic = "+proj=longlat +datum=WGS84 +no_defs +type=crs";
const cesiumCartesian =
  "+proj=geocent +datum=WGS84 +units=m +no_defs +type=crs";

const cesiumCartesianToCartographicTransform = proj4(
  cesiumCartesian,
  cesiumCartographic
);

export function getNewellsNormal(points: Vector3Like[]) {
  // find normal with Newell's method
  let n = [0.0, 0.0, 0.0];

  for (let i = 0; i < points.length; i++) {
    let nex = i + 1;

    if (nex == points.length) {
      nex = 0;
    }

    n[0] = n[0] + (points[i].y - points[nex].y) * (points[i].z + points[nex].z);
    n[1] = n[1] + (points[i].z - points[nex].z) * (points[i].x + points[nex].x);
    n[2] = n[2] + (points[i].x - points[nex].x) * (points[i].y + points[nex].y);
  }

  let b = new Vector3(n[0], n[1], n[2]);

  return b.normalize();
}

export function projectPolygonToPlane(
  polygon3D: [number, number, number][],
  polygonHoles3D: [number, number, number][]
) {
  const points: Vector3[] = [];
  for (const point3D of polygon3D) {
    points.push(new Vector3(...point3D));
  }

  const holePoints: Vector3[] = [];
  for (const holePoint3D of polygonHoles3D) {
    holePoints.push(new Vector3(...holePoint3D));
  }

  const allPoints = points.concat(holePoints);

  const normal = getNewellsNormal(allPoints);

  // 4. Create the quaternion to rotate the plane so that its normal aligns with (0, 0, 1).
  const quaternion = new Quaternion().setFromUnitVectors(
    normal,
    new Vector3(0, 0, 1)
  );

  const matrix = new Matrix4().makeRotationFromQuaternion(quaternion);

  // 5. Project the points using the quaternion.
  const projected = points.map((point) => point.clone().applyMatrix4(matrix));

  return {
    projected,
    matrix,
  };
}

export function triangulate3DPolygon(rings: [number, number, number][][]) {
  if (!rings || rings.length === 0) return null;

  // maybe even scale down to range of [-1, 1]
  const centralizedRings = rings;

  let outerRing = centralizedRings[0];
  if (outerRing.length === 0) return null;

  // Compute the projection for the outer ring.
  const projection = projectPolygonToPlane(
    outerRing,
    centralizedRings.slice(1).flat()
  );

  const outerPolygons = projection.projected.flatMap((p) => [p.x, p.y]);
  const holes: number[] = [];

  let lastHoleIndex = centralizedRings[0].length;

  const holeIndices: number[] = [];

  for (let i = 1; i < centralizedRings.length; i++) {
    let hole = centralizedRings[i];
    if (hole.length === 0) continue;
    holeIndices.push(lastHoleIndex);
    lastHoleIndex += hole.length;
    const hole3d = hole.map((pt) =>
      pt.length < 3
        ? new Vector3(pt[0], pt[1], 0)
        : new Vector3(pt[0], pt[1], pt[2])
    );
    hole3d.forEach((pt) => {
      pt.applyMatrix4(projection.matrix);
      holes.push(pt.x, pt.y);
    });
  }

  const points2d = outerPolygons.concat(holes);

  const triangles = earcut.default(points2d, holeIndices, 2);

  return triangles;
}

export function getBBoxesFromMeshes(
  meshes: { position: Float32Array }[],
  globalTransformPoint: Vector3
) {
  const cartesianBox = new Box3();
  const cartographicBox = new Box3();

  for (const mesh of meshes) {
    for (let i = 0; i < mesh.position.length; ) {
      const [x, y, z] = [
        mesh.position[i++] + globalTransformPoint.x,
        mesh.position[i++] + globalTransformPoint.z,
        mesh.position[i++] - globalTransformPoint.y,
      ];
      const cartesian = [x, -z, y]; // back to not y up coordinate system
      const [longitude, latitude, height] =
        cesiumCartesianToCartographicTransform.forward(cartesian);

      cartesianBox.expandByPoint(new Vector3(x, y, z));
      cartographicBox.expandByPoint(
        new Vector3(
          (longitude / 180) * Math.PI,
          (latitude / 180) * Math.PI,
          height
        )
      );
    }
  }

  return { cartesianBox, cartographicBox };
}

export function convertEPSGFromCityJSONToProj4(potentialEPSG?: string) {
  if (!potentialEPSG)
    throw new Error(
      "Source coordinate system is undefined, please provice an EPSG code!"
    );

  let epsg: string | undefined = undefined;

  try {
    epsg = potentialEPSG
      .split(",")
      .filter((part) => part.includes("EPSG"))
      .map((part) => part.split(":").pop())[0];
  } catch {}
  try {
    epsg = potentialEPSG.split("/").slice(-1)[0];
  } catch {}

  if (!epsg)
    throw new Error("Could not parse EPSG code please provide a valid one!");

  if (!proj4List[`EPSG:${epsg}`]?.[1])
    throw new Error(`Did not find "${epsg}" in epsg list!`);

  return proj4List[`EPSG:${epsg}`][1];
}
