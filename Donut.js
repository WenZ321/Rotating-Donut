(function () {
  // Get the canvas and its 2D context
  const canvas = document.getElementById("donutCanvas");
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  // ======================
  // DONUT
  // ======================
  // R: distance from donut center to the center of the tube
  // r: radius of the tube
  // stepsPhi: subdivisions around the tube itself
  // stepsTheta: subdivisions around the central hole
  const R = 150;       
  const r = 40;        
  const stepsPhi = 20; 
  const stepsTheta = 20;

  // The focal length is used for perspective projection:
  // A larger focal length leads to a weaker perspective effect (things don't shrink as much in the distance),
  // while a smaller focal length exaggerates perspective.
  const focalLength = 1000;

  // =================================
  // Generating Donut
  // =================================
  // We'll store them in a 2D array: donutPoints[phiIndex][thetaIndex] = [x, y, z]
  // Parametric equations for a torus:
  //   x = (R + r*cos(phi))*cos(theta)
  //   y = (R + r*cos(phi))*sin(theta)
  //   z = r*sin(phi)
  // where phi goes around the small tube, and theta goes around the central hole.
  let donutPoints = new Array(stepsPhi);
  for (let i = 0; i < stepsPhi; i++) {
    // Create the inner array for each phi step
    donutPoints[i] = new Array(stepsTheta);
    
    // phi is the angle for the circular cross section of the tube
    const phi = (i / stepsPhi) * 2 * Math.PI;
    const cosPhi = Math.cos(phi);
    const sinPhi = Math.sin(phi);

    for (let j = 0; j < stepsTheta; j++) {
      // theta is the angle around the donut's central hole
      const theta = (j / stepsTheta) * 2 * Math.PI;
      const cosTheta = Math.cos(theta);
      const sinTheta = Math.sin(theta);

      // Calculate 3D coordinates using the donut parametric formula
      const x = (R + r * cosPhi) * cosTheta;
      const y = (R + r * cosPhi) * sinTheta;
      const z = r * sinPhi;

      // Store the 3D point
      donutPoints[i][j] = [x, y, z];
    }
  }

  // Track the animation angle in radians.
  // We'll rotate the donut by this angle around X and Y axes.
  let angle = 0;

  // ======================
  // MAIN ANIMATION LOOP
  // ======================
  function animate() {
    angle += 0.01;
    draw(angle);
    requestAnimationFrame(animate);
  }

  /**
   * Draw the donut at a given rotation angle.
   * @param {number} angle - Current rotation angle in radians.
   */
  function draw(angle) {
    // Clear the entire canvas so we can draw the new frame
    ctx.clearRect(0, 0, width, height);

    // ======================
    // ROTATION 
    // ======================
    // Composition of rotating around Y and then around X
    // We'll apply it to every donut vertex.
    const rotationMatrix = multiplyMatrix(rotationY(angle),rotationX(angle));

    // Array to hold the 2D coordinates of the projected points after they've been rotated and projected.
    let projected2D = new Array(stepsPhi);

    // For each (phi, theta) point:
    for (let i = 0; i < stepsPhi; i++) {
      projected2D[i] = new Array(stepsTheta);
      
      for (let j = 0; j < stepsTheta; j++) {

        const pt3D = donutPoints[i][j];

        // Rotate the point using our rotation matrix
        const rotated = multiplyMatrixVector(rotationMatrix, pt3D);

        // ======================
        // PERSPECTIVE PROJECTION
        // ======================
        // We treat the focal length as if the "camera" is located at z = -focalLength.
        // The further a point is from the camera, the smaller it appears.
        const zOffset = rotated[2] + focalLength; 
        const scale = focalLength / zOffset;

        // Convert 3D -> 2D using perspective scaling,
        // then shift the coordinate so (0,0) is in the middle of the canvas.
        const x2d = rotated[0] * scale + width / 2;
        const y2d = rotated[1] * scale + height / 2;

        // Store 2D coordinate for later drawing
        projected2D[i][j] = { x: x2d, y: y2d };
      }
    }

    // ================================
    // DRAWING DONUTT
    // ================================
    ctx.strokeStyle = "#FF0000"; // Color of donut 
    
    // For each point, connect it to its neighbors in phi and theta directions to form a grid
    for (let i = 0; i < stepsPhi; i++) {
      for (let j = 0; j < stepsTheta; j++) {
        // Wrap around using modulo so the last point connects back to the first 
        const nextI = (i + 1) % stepsPhi;
        const nextJ = (j + 1) % stepsTheta;

        // Current point in 2D
        const p = projected2D[i][j];
        // The neighbor in the "phi" direction
        const pPhi = projected2D[nextI][j];
        // The neighbor in the "theta" direction
        const pTheta = projected2D[i][nextJ];

        // Draw line to the neighbor in the phi direction
        ctx.beginPath();
        ctx.moveTo(p.x, p.y);
        ctx.lineTo(pPhi.x, pPhi.y);
        ctx.stroke();

        // Draw line to the neighbor in the theta direction
        ctx.beginPath();
        ctx.moveTo(p.x, p.y);
        ctx.lineTo(pTheta.x, pTheta.y);
        ctx.stroke();
      }
    }
  }

  /**
   * Rotation around the X-axis
   *
   * R_x(ax) =
   * [  1    0       0     ]
   * [  0   cos(ax) -sin(ax)]
   * [  0   sin(ax)  cos(ax)]
   */
  function rotationX(ax) {
    const c = Math.cos(ax);
    const s = Math.sin(ax);
    return [
      [1,  0,  0],
      [0,  c, -s],
      [0,  s,  c]
    ];
  }

  /**
   *Rotation around the Y-axis 
   *
   * R_y(ay) =
   * [  cos(ay)  0   sin(ay) ]
   * [     0     1     0    ]
   * [ -sin(ay)  0   cos(ay) ]
   */
  function rotationY(ay) {
    const c = Math.cos(ay);
    const s = Math.sin(ay);
    return [
      [ c,  0,  s],
      [ 0,  1,  0],
      [-s,  0,  c]
    ];
  }

  function multiplyMatrix(A, B) {
    let result = [];
    for (let i = 0; i < 3; i++) {
      result[i] = [];
      for (let j = 0; j < 3; j++) {
        let sum = 0;
        for (let k = 0; k < 3; k++) {
          sum += A[i][k] * B[k][j];
        }
        result[i][j] = sum;
      }
    }
    return result;
  }


  function multiplyMatrixVector(M, v) {
    let x = M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2];
    let y = M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2];
    let z = M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2];
    return [x, y, z];
  }


  animate();
})();
