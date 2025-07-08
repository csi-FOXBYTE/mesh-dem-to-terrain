// Sutherland-Hodgeman polygon clipping algorithm

export function polygonclip(
  points: [number, number, number][],
  bbox: [number, number, number, number]
) {
  let prev: [number, number, number],
    prevInside: boolean,
    i: number,
    p: [number, number, number],
    inside: boolean,
    result: [number, number, number][] = [];

  // clip against each side of the clip rectangle
  for (let edge = 1; edge <= 8; edge *= 2) {
    result = [];
    prev = points[points.length - 1];
    prevInside = !(bitCode(prev, bbox) & edge);

    for (i = 0; i < points.length; i++) {
      p = points[i];
      inside = !(bitCode(p, bbox) & edge);

      // if segment goes through the clip window, add an intersection
      if (inside !== prevInside) result.push(intersect(prev, p, edge, bbox));

      if (inside) result.push(p); // add a point if it's inside

      prev = p;
      prevInside = inside;
    }

    points = result;

    if (!points.length) break;
  }

  return result;
}

// intersect a segment against one of the 4 lines that make up the bbox

function intersect(
  a: [number, number, number],
  b: [number, number, number],
  edge: number,
  bbox: [number, number, number, number]
) {
  let t: number, x: number, y: number, z: number;

  if (edge & 8) {
    if (b[1] === a[1]) {
      throw new Error("Segment is horizontal; cannot intersect with top clip edge.");
    }
    // top: y = bbox[3]
    t = (bbox[3] - a[1]) / (b[1] - a[1]);
    x = a[0] + t * (b[0] - a[0]);
    y = bbox[3];
  } else if (edge & 4) {
    if (b[1] === a[1]) {
      throw new Error("Segment is horizontal; cannot intersect with bottom clip edge.");
    }
    // bottom: y = bbox[1]
    t = (bbox[1] - a[1]) / (b[1] - a[1]);
    x = a[0] + t * (b[0] - a[0]);
    y = bbox[1];
  } else if (edge & 2) {
    if (b[0] === a[0]) {
      throw new Error("Segment is vertical; cannot intersect with right clip edge.");
    }
    // right: x = bbox[2]
    t = (bbox[2] - a[0]) / (b[0] - a[0]);
    x = bbox[2];
    y = a[1] + t * (b[1] - a[1]);
  } else if (edge & 1) {
    if (b[0] === a[0]) {
      throw new Error("Segment is vertical; cannot intersect with left clip edge.");
    }
    // left: x = bbox[0]
    t = (bbox[0] - a[0]) / (b[0] - a[0]);
    x = bbox[0];
    y = a[1] + t * (b[1] - a[1]);
  } else {
    throw new Error("Invalid edge value");
  }

  if (t > 1 || t < 0) throw new Error("Calculated parameter t is out of range.")

  z = a[2] + t * (b[2] - a[2]);

  return [x, y, z] as [number, number, number];
}

// bit code reflects the point position relative to the bbox:

//         left  mid  right
//    top  1001  1000  1010
//    mid  0001  0000  0010
// bottom  0101  0100  0110

function bitCode(
  p: [number, number, number],
  bbox: [number, number, number, number]
) {
  let code = 0;

  if (p[0] < bbox[0]) code |= 1; // left
  else if (p[0] > bbox[2]) code |= 2; // right

  if (p[1] < bbox[1]) code |= 4; // bottom
  else if (p[1] > bbox[3]) code |= 8; // top

  return code;
}
