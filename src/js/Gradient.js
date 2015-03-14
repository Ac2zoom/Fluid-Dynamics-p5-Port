/**
 * Created by Akshay on 3/13/2015.
 */
var Gradient = function(){
    this.DEFAULT_BLACK_TO_WHITE = 0;
    this.DEFAULT_RANDOM = 1;
    this.DEFAULT_SPECTRUM = 2;
    this.DEFAULT_INFRARED = 3;
    this.DEFAULT_BLACKBODY = 4;
    this.DEFAULT_NEON = 5;
    this.DEFAULT_WINTER = 6;
    this.DEFAULT_SUMMER = 7;

    var GradientNode = function(l, c) {
        this.location = l;
        this.col = c;
    };

    var nodes = [];

    this.addNode = function(location, col) {
        if ( nodes.length == 0 ) {
            nodes[0] = new GradientNode(location, col);
        }
        else {
            this.newNodes = nodes;
            this.newNodes[nodes.length] = new GradientNode(location, col);
            this.newNodes[nodes.length+1] = new GradientNode();
            nodes = this.newNodes;
        }
    };

    this.clear = function() {
        nodes = [];
    };

    this.makeDefaultGradient = function(defaultGradient) {

        this.clear();

        switch ( defaultGradient ) {
            case this.DEFAULT_BLACK_TO_WHITE:
                this.addNode ( 0, 0x000000 );
                this.addNode ( 1, 0xFFFFFF );
                break;
            case this.DEFAULT_RANDOM:
                this.makeRandomGradient(4);
                break;
            case this.DEFAULT_SPECTRUM:
                this.addNode ( 0,    0xFF0000 );
                this.addNode ( 0.25, 0xFFFF00 );
                this.addNode ( 0.5,  0x00FF00 );
                this.addNode ( 0.75, 0x00FFFF );
                this.addNode ( 1,    0x0000FF );
                break;
            case this.DEFAULT_INFRARED:
                this.addNode ( 0,    0x000000 );
                this.addNode ( 1.0/6,0x000080 );
                this.addNode ( 2.0/6,0x800080 );
                this.addNode ( 3.0/6,0xFF0000 );
                this.addNode ( 4.0/6,0xFF8000 );
                this.addNode ( 5.0/6,0xFFFF00 );
                this.addNode ( 1,    0xFFFFFF );
                break;
            case this.DEFAULT_BLACKBODY:
                this.addNode ( 0,    0x000000 );
                this.addNode ( 1.0/5,0x0040FF );
                this.addNode ( 2.0/5,0x00C0FF );
                this.addNode ( 3.0/5,0xFF4000 );
                this.addNode ( 4.0/5,0xFFC000 );
                this.addNode ( 1,    0xFFFFFF );
                break;
            case this.DEFAULT_NEON:
                this.addNode ( 0,    0x000000 );
                this.addNode ( 0.25, 0x3333FF );
                this.addNode ( 0.5,  0x0099FF );
                this.addNode ( 0.75, 0xE60080 );
                this.addNode ( 1,    0xFF00FF );
                break;
            case this.DEFAULT_WINTER:
                this.addNode ( 0,    0x4C80FF );
                this.addNode ( 0.5,  0xE6E6E6 );
                this.addNode ( 1,    0x999999 );
                break;
            case this.DEFAULT_SUMMER:
                this.addNode ( 0,    0x334CFF );
                this.addNode ( 0.25, 0xFF0080 );
                this.addNode ( 0.5,  0xFF8033 );
                this.addNode ( 0.75, 0xCC4C00 );
                this.addNode ( 1,    0xFFCC00 );
                break;
        }
    };

    this.makeRandomGradient = function(numColours) {
        var location = Number();
        var locationMin = Number();
        var locationMax = Number();
        var r = Number();
        var g = Number();
        var b = Number();

        this.clear();

        for(var n = 0; n < numColours; n++) {
            if ( n == 0 ) {
                location = 0.0;
            } else if ( n == numColours-1 ){
                location = 1.0;
            } else {
                locationMin = float(n)/numColours;
                locationMax = float(n+1)/numColours;
                location = (Math.random() * locationMax) + locationMin;
            }

            r = parseInt(Math.random() * 2.5) * 128;
            g = parseInt(Math.random() * 2.5) * 128;
            b = parseInt(Math.random() * 2.5) * 128;

            this.addNode(location, color(r,g,b));
        }
    };

    this.makeDefaultGradient(this.DEFAULT_BLACK_TO_WHITE);

    this.getColour = function(location) {
        var bandLocation = Number(),
            bandScale = Number(),
            bandDelta = Number(),
            r = Number(),
            g = Number(),
            b = Number();

        for ( var c = 0; c < nodes.length-1; c++ ) {
            if ( location >= nodes[c].location && location <= nodes[c+1].location ) {
                bandScale = nodes[c+1].location - nodes[c].location;
                bandLocation = location - nodes[c].location;
                bandDelta = bandLocation / bandScale;

                r = bandDelta * ( red(nodes[c+1].col) - red(nodes[c].col) ) + red(nodes[c].col);
                g = bandDelta * ( green(nodes[c+1].col) - green(nodes[c].col) ) + green(nodes[c].col);
                b = bandDelta * ( blue(nodes[c+1].col) - blue(nodes[c].col) ) + blue(nodes[c].col);
                return color(r,g,b);
            }
        }
        return color(0,0,0);
    };

    this.makeArrayOfColours = function(numColours) {

        var location = Number(),
            bandLocation = Number(),
            bandScale = Number(),
            bandDelta = Number(),
            cols = [];

        for ( var i = 0; i < numColours; i++ ) {
            location = float(i) / (numColours-1);
            cols[i] = this.getColour(location);
        }

        return cols;
    }
};
