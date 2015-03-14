var n = 128;
var gridSize = n+2;
var pixelSize = 2;

var u = [[gridSize],[gridSize]];
var v = [[gridSize],[gridSize]];
var uPrev = [[gridSize],[gridSize]];
var vPrev = [[gridSize],[gridSize]];
var dens = [[gridSize],[gridSize]];
var densPrev = [[gridSize],[gridSize]];

var viscosity = 0.0001;
var dt = 0.2;
var diff = 0.0001;

var prevMouseX;
var prevMouseY;

var showDensity = true;
var showVelocity = false;
var showNormals = false;
var showLighting = true;
var showGradient = false;

var shiftPressed = false;

var gradient = new Gradient();
var gradientColours = [];
var numGradientColours = 256;

var densityBrushSize = n/10;
var velocityBrushSize = n/20;  // Ditto velocity
var lineSpacing = n/20;

function setup() {
    createCanvas(n*pixelSize, n*pixelSize);
    background(0);
    noStroke();

    gradient.makeDefaultGradient(gradient.DEFAULT_BLACKBODY);
    gradientColours = gradient.makeArrayOfColours(numGradientColours);

    reset();

    frameRate(1000);
}

function draw() {
    setForce();
    calculateVelocity(u, v, uPrev, vPrev, viscosity, dt);
    calculateDensity(dens, densPrev, u, v, diff, dt);
    drawField(dens);
}

var reset = function() {
    initVelocity();
    initDensity();
};

var initField = function(f) {
    for ( var i = 0; i < gridSize; i++ ) {
        for ( var j = 0; j < gridSize; j++ ) {
            f[i][j] = 0.0;
        }
    }
};

var initVelocity = function() {
    initField(u);
    initField(v);
    initField(uPrev);
    initField(vPrev);
};

var initDensity = function() {
    initField(dens);
    initField(densPrev);
};

var ix = function(i, j) {
    return i + j*(n+2);
};

var addSource = function(x, s, dt) {
    for ( var i = 0; i < gridSize; i++ ) {
        for ( var j = 0; j < gridSize; j++ ) {
            x[i][j] += s[i][j]*dt;
        }
    }
};

var setBnd = function(b, x) {

    var i;
    for ( i = 1; i <= n; i++ ) {
        if ( b==1 )  x [0  ][i  ] = -x [1  ][i  ];  else  x [0  ][i  ] = x [1  ][i  ];
        if ( b==1 )  x [n+1][i  ] = -x [n  ][i  ];  else  x [n+1][i  ] = x [n  ][i  ];
        if ( b==2 )  x [i  ][0  ] = -x [i  ][1  ];  else  x [i  ][0  ] = x [i  ][1  ];
        if ( b==2 )  x [i  ][n+1] = -x [i  ][n  ];  else  x [i  ][n+1] = x [i  ][n  ];
    }
    x [0  ][0  ] = 0.5 * ( x [1  ][0  ] + x [0  ][1  ] );
    x [0  ][n+1] = 0.5 * ( x [1  ][n+1] + x [0  ][n  ] );
    x [n+1][0  ] = 0.5 * ( x [n  ][0  ] + x [n+1][1  ] );
    x [n+1][n+1] = 0.5 * ( x [n  ][n+1] + x [n+1][n  ] );
};

var diffuse = function (b, x, x0, diff, dt) {
    var i = Number(),
        j = Number(),
        k = Number(),
        a = dt * diff * n * n;
    for ( k = 0; k < 20; k++ ) {
        for ( i = 1; i <= n; i++ ) {
            for ( j = 1; j <= n; j++ ) {
                x[i][j] = (x0[i][j] + a * (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1])) / (1+4*a);
            }
        }
        setBnd(b, x);
    }
};

var project = function(u, v, p, div) {
    var i = Number(),
        j = Number(),
        k = Number(),
        h = Number();

    h = 1.0/n;
    for ( i = 1; i<=n; i++ ) {
        for ( j = 1; j <= n; j++ ) {
            div [i][j] = -0.5 * h * ( u [i+1][j] - u [i-1][j] + v [i][j+1] - v [i][j-1] );
            p [i][j] = 0;
        }
    }
    setBnd( 0, div );
    setBnd( 0, p );

    for ( k=0; k<20; k++ ) {
        for ( i=1; i<=n; i++ ) {
            for ( j=1; j<=n; j++ ) {
                p [i][j] = ( div [i][j] + p [i-1][j] + p [i+1][j] + p [i][j-1] + p [i][j+1] ) / 4;
            }
        }
        setBnd ( 0, p );
    }

    for ( i=1; i<=n; i++ ) {
        for ( j=1; j<=n; j++ ) {
            u [i][j] -= 0.5 * ( p [i+1][j] - p [i-1][j] ) / h;
            v [i][j] -= 0.5 * ( p [i][j+1] - p [i][j-1] ) / h;
        }
    }
    setBnd ( 1, u );
    setBnd ( 2, v );
};

var advect = function( b, d, d0, u, v, dt) {
    var i = Number(),
        j = Number(),
        i0 = Number(),
        j0 = Number(),
        i1 = Number(),
        j1 = Number(),
        x = Number(),
        y = Number(),
        s0 = Number(),
        t0 = Number(),
        s1 = Number(),
        t1 = Number(),
        dt0 = Number();

    dt0 = dt*n;
    for ( i=1; i<=n; i++ ) {
        for ( j=1; j<=n; j++ ) {
            x = i - dt0 * u [i][j];
            y = j - dt0 * v [i][j];

            x = Math.max(0.5, x);
            x = Math.min(n + 0.5, x);

            //i0 = (int)x;
            i0 = Math.floor(x);
            i1 = i0 + 1;

            y = Math.max(0.5, y);
            y = Math.min(n + 0.5, y);

            //j0 = (int)j;
            j0 = Math.floor(y);
            j1 = j0 + 1;

            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;

            d [i][j] = s0 * ( t0 * d0 [i0][j0] + t1 * d0 [i0][j1] ) +
            s1 * ( t0 * d0 [i1][j0] + t1 * d0 [i1][j1] );
        }
    }
    setBnd(b, d);
};

var calculateVelocity = function(u, v, u0, v0, visc, dt) {
    addSource ( u, u0, dt );
    addSource ( v, v0, dt );

    var tmp;
    tmp = u;  u = u0;  u0 = tmp;
    tmp = v;  v = v0;  v0 = tmp;

    diffuse ( 1, u, u0, visc, dt );
    diffuse ( 2, v, v0, visc, dt );

    project ( u, v, u0, v0 );

    tmp = u;  u = u0;  u0 = tmp;
    tmp = v;  v = v0;  v0 = tmp;

    advect ( 1, u, u0, u0, v0, dt );
    advect ( 2, v, v0, u0, v0, dt );

    project ( u, v, u0, v0 );
};

var calculateDensity = function(x, x0, u, v, diff, dt) {
    var tmp;

    addSource ( x, x0, dt );
    tmp = x; x = x0; x0 = tmp;
    diffuse ( 0, x, x0, diff, dt );
    tmp = x; x = x0; x0 = tmp;
    advect ( 0, x, x0, u, v, dt );
};

var drawField = function(f) {
    var x = Number(),
        y = Number(),
        s = Number(),
        col,
        ax = Number(),
        ay = Number(),
        az = Number(),
        bx = Number(),
        by = Number(),
        bz = Number(),
        nx = 0,
        ny = 0,
        nz = 0,
        vu = Number(),
        vv = Number(),
        l = Number(),
        lx = 0,
        ly = 0,
        lz = 0,
        vx = 0,
        vy = 0,
        vz = 0,
        hx = 0,
        hy = 0,
        hz = 0,
        ndoth = Number(),
        ndot1 = Number(),
        w = Number(),
        d = Number();

    background(0);

    for (y = 1; y <= n; y++ ) {
        for ( x = 1; x <= n; x++ ) {

            d = dens[x][y];

            if ( showNormals == true || showLighting == true) {
                //ax = 1;
                //ay = 0;
                az = d - dens[x+1][y];
                //bx = 0;
                //by = 1;
                bz = d - dens[x][y+1];

                //nx = ay*bz - az*by;
                //ny = az*bx - ax*bz;
                //nz = ax*by - ay*bx;

                nx = -az;
                ny = -bz;
                nz = 1;

                //l = sqrt(nx*nx + ny*ny + nz*nz); if ( l != 0 ) { nx /= l; ny /= l; nz /= l; }
                l = -sqrt(nx*nx + ny*ny + 1); nx /= l; ny /= l; nz /= l;
            }

            if ( showDensity == true) {
                s = abs(range(numGradientColours * d,0,numGradientColours-1));
                noStroke();

                col = gradientColours[parseInt(s)];

                if ( showLighting == true) {
                    /*
                     lx = n/2 - x;
                     ly = n/2 - y;
                     lz = 1000 - 1000*dens[x,y];
                     l = sqrt(lx*lx + ly*ly + lz*lz); if ( l != 0 ) { lx /= l; ly /= l; lz /= l; }

                     vx = n/2 - x;
                     vy = n/2 - y;
                     vz = n;
                     l = sqrt(vx*vx + vy*vy + vz*vz); if ( l != 0 ) { vx /= l; vy /= l; vz /= l; }

                     lx = 0; ly = 0; lz = 1;
                     vx = 0; vy = 0; vz = 1;

                     hx = (lx + vx)/2;
                     hy = (ly + vy)/2;
                     hz = (lz + vz)/2;
                     l = sqrt(hx*hx + hy*hy + hz*hz); if ( l != 0 ) { hx /= l; hy /= l; hz /= l; }

                     // ndoth = nx*hx + ny*hy + nz*hz;
                     */

                    ndoth = nz;

                    w = ((s*512)/numGradientColours) * pow(ndoth, 100000);

                    col = blendColor ( col, color(w), ADD );

                }

                fill ( col );
                rect ( (x-1) * pixelSize, (y-1) * pixelSize, pixelSize, pixelSize );
            }

            if ( showVelocity == true) {
                if ( (x % lineSpacing) == 0 && (y % lineSpacing) == 0 ) {
                    noFill();
                    stroke(255,0,0);
                    vu = range(500 * u[x][y], -50,50);
                    vv = range(500 * v[x][y], -50,50);
                    line( (x-1)*pixelSize, (y-1)*pixelSize, (x-1)*pixelSize + vu, (y-1)*pixelSize + vv);
                }
            }

            if ( showNormals == true) {
                if ( (x % lineSpacing) == 0 && (y % lineSpacing) == 0 ) {
                    noFill();
                    stroke(255);

                    vu = range(500 * nx, -50,50);
                    vv = range(500 * ny, -50,50);
                    line( (x-1)*pixelSize, (y-1)*pixelSize, (x-1)*pixelSize + vu, (y-1)*pixelSize + vv);
                }
            }
        }
    }

    if ( showGradient == true) {
        noStroke();
        for ( var i=0; i < numGradientColours; i++ ) {
            fill(gradientColours[i]);
            rect ( i % width, i / width, 1,10 );
        }
    }
};

var setForce = function() {
    initField(densPrev);
    initField(uPrev);
    initField(vPrev);

    if ( mousePressed == true) {
        var x = Number(),
            y = Number();

        x = (mouseX * n) / (width) + 1;
        y = (mouseY * n) / (height) + 1;

        if (keyPressed) {
            shiftPressed = false;
        }

        if ( mouseButton == LEFT && shiftPressed !== true) {
            setForceArea ( uPrev, x, y, mouseX - prevMouseX, velocityBrushSize );
            setForceArea ( vPrev, x, y, mouseY - prevMouseY, velocityBrushSize );
        } else {
            var m = (mouseX - prevMouseX) + (mouseY - prevMouseY);
            setForceArea ( densPrev, x, y, range(abs(m),0,2), densityBrushSize );
        }
    }
    prevMouseX = mouseX;
    prevMouseY = mouseY;
};

var range = function( f, minf, maxf) {
    return Math.max(Math.min(f, maxf), minf);
};

var setForceArea = function(field, x, y, s, r) {
        var i = Number(),
            j = Number(),
            dx = Number(),
            dy = Number(),
            f = Number();

        for ( i = parseInt(range(x-r,1,n)); i <= parseInt(range(x+r, 1, n)); i++) {
            dx = x - i;
            for (j = parseInt(range(y-r, 1, n)); j <= parseInt(range(y+r, 1, n)); j++) {
                dy = y - j;
                f = 1 - ( sqrt(dx*dx + dy*dy) / r);
                field[i][j] += range(f,0,1) * s;
            }
    }
};

function keyPressed() {
    if ( keyCode == 32 ) {
        gradient.makeRandomGradient(4);
        gradientColours = gradient.makeArrayOfColours(numGradientColours);
    }

    if ( key == 'v' ){
        showVelocity = !showVelocity;
    }
    if ( key == 'n' ){
        showNormals = !showNormals;
    }
    if ( key == 'd' ) {
        showDensity = !showDensity;
    }
    if ( key == 'l' ){
        showLighting = !showLighting;
    }
    if ( key == 'r' ){
        reset();
    }
    if ( key == '+' ){
        viscosity += 0.0001;
    }
    if ( key == '-' ){
        viscosity -= 0.0001;
    }
    viscosity = Math.max(viscosity, 0);

    if ( key >= '0' && key <= '7' ) {
        gradient.makeDefaultGradient(key - '0');
        gradientColours = gradient.makeArrayOfColours(numGradientColours);
    }

    if ( key == 'g' ) {
        showGradient = ! showGradient;
    }
    if ( key == 's' ) {
        save("screenshot-" + year() + month() + day() + hour() + minute() + second() + ".png");
    }
    if ( keyCode == SHIFT ) {
        shiftPressed = true;
    }
};
