/*
 * L.Transformation is an utility class to perform simple point transformations through a 2d-matrix.
 */


L.Matrix = function () {
	this.padding = -1;
    this.m = this._newM();
};

L.Matrix.prototype = {
	_newM : function(){
        return [[this.padding],
            [this.padding, 1,0,0],
            [this.padding, 0,1,0],
            [this.padding, 0,0,1]];
    },

    _cloneM : function(m){
        var clone = this._newM();
        for (var x = 1; x<4; ++x)
        {
                for (var y = 1; y<4; ++y)
                {
                        clone[x][y] = m[x][y];
                }
        }
        return clone;
    },

    /**
    * Multiplies matrix with given matrix and returns resulting matrix
    *
    * @param o {L.Matrix}
    * @returns {L.Matrix}
    */
    mul : function(o) {
        var product = new L.Matrix();

        for (var x = 1; x<4; ++x)
        {
            for (var y = 1; y<4; ++y)
            {
                var sum = 0;
                for (var z = 1; z<4; ++z){
                    sum += this.m[x][z] * o.m[z][y];
                }
                product.m[x][y] = sum;
            }
        }
        this.m = product.m;
        return this;
    },

    determinant : function() {
        var min = this.minors(this.m);
        var a = this.m[1][1] * min[1][1];
        var b = this.m[1][2] * min[1][2];
        var c = this.m[1][3] * min[1][3];
        return (a - b + c);
    },
    
    determinant2 : function(a, b , c, d) {
        return (a * d) - (b * c);
    },
    
    minors : function(m) {
        var minor11 = this.determinant2(m[2][2], m[2][3], m[3][2], m[3][3]) ;
        var minor12 = this.determinant2(m[2][1], m[2][3], m[3][1], m[3][3]) ;
        var minor13 = this.determinant2(m[2][1], m[2][2], m[3][1], m[3][2]) ;
        var minor21 = this.determinant2(m[1][2], m[1][3], m[3][2], m[3][3]) ;
        var minor22 = this.determinant2(m[1][1], m[1][3], m[3][1], m[3][3]) ;
        var minor23 = this.determinant2(m[1][1], m[1][2], m[3][1], m[3][2]) ;
        var minor31 = this.determinant2(m[1][2], m[1][3], m[2][2], m[2][3]) ;
        var minor32 = this.determinant2(m[1][1], m[1][3], m[2][1], m[2][3]) ;
        var minor33 = this.determinant2(m[1][1], m[1][2], m[2][1], m[2][3]) ; 
        var padding = -1;
        return  [[padding],
                [padding, minor11,minor12,minor13],
                [padding, minor21,minor22,minor23],
                [padding, minor31,minor32,minor33]];
    },
    
    cofactors : function(m) {
        var c = this._cloneM(m);
        c[1][2] = m[1][2] * -1;
        c[2][1] = m[2][1] * -1;
        c[2][3] = m[2][3] * -1;
        c[3][2] = m[3][2] * -1;
        return c;
    },
    
    adjugate: function(m) {
        var a = this._cloneM(m);
        a[1][2] = m[2][1];
        a[2][1] = m[1][2];
        a[1][3] = m[3][1];
        a[3][1] = m[1][3];
        a[2][3] = m[3][2];
        a[3][2] = m[2][3];
        return a;
    },
    
    inverse : function() {
        var det = this.determinant();
        // we assume that the matrix is always invertible, so wrong but restful :)
        // well, it would throw a division by 0 exception a bit later, thats it
        var m = this.adjugate(this.cofactors(this.minors(this.m)));
        
        var inverse = new L.Matrix();
        for (var x = 1; x<4; ++x)
        {
                for (var y = 1; y<4; ++y)
                {
                        inverse.m[x][y] = m[x][y] * (1/det);
                }
        }
        return inverse;
    }
};

L.Transform =  function() {
    this.m = new L.Matrix(); 
};

L.Transform.prototype = {

    inverse: function() {
        var inverseM = this.m.inverse();
        var inverse = new L.Transform();
        inverse.m = inverseM;
        return inverse;
    },
    
    multiply: function(t){
        if(t instanceof L.Matrix){
            this.m.mul(t);
        }
        else{
            this.m.mul(t.m);
        }
        return this;
    },

    translate: function(tx,ty) {
        var transMat = new L.Matrix();
        transMat.m[3][1] = tx;
        transMat.m[3][2] = ty;
        this.m.mul(transMat);
        return this;
    },

    /**
    * Returns the current translate of the transformation matrix
    *
    * @returns {L.Point}
    */
    getTranslate: function(){
        return new L.Point(this.m.m[3][1], this.m.m[3][2]);
    },

    /**
    * Resets the translation of the matrix to the given x and y
    * or 0, 0 if no points were provided. Returns itself
    * 
    * @param x int
    * @param y int
    * @returns {L.Transform}
    */
    resetTranslate: function(x,y) {            
        this.m.m[3][1] = x || 0;
        this.m.m[3][2] = y || 0;
        return this;
    },

    /**
    * Scales with given scale on x-axis and 
    * given scale on y-axis, around given origin
    * 
    * If no sy is provided the scale will be proportional
    * @param sx Number
    * @param sy Number
    * @param origin {L.Point}|{}
    * @returns {L.Transform}
    */
    scale: function(sx, sy, origin) {
        var scaleMat = new L.Matrix();
        if(origin !== undefined)
        {
            var tr1 = new L.Matrix();
            tr1.m[3][1] = -origin.x;
            tr1.m[3][2] = -origin.y;
            scaleMat.mul(tr1);

            var tr2 = new L.Matrix();
            tr2.m[1][1] = sx;
            tr2.m[2][2] = sy;
            scaleMat.mul(tr2);

            var tr3 = new L.Matrix();
            tr3.m[3][1] = origin.x;
            tr3.m[3][2] = origin.y;
            scaleMat.mul(tr3);
        }
        else
        {
            scaleMat.m[1][1] = sx;
            scaleMat.m[2][2] = sy;
        }
        this.m.mul(scaleMat);
        return this;
    },

    rotate: function(r, origin){
        var rotMat = new L.Matrix();
        var rGrad = r * Math.PI / 180.0;
        var cosR = Math.cos(rGrad);
        var sinR = Math.sin(rGrad);

        if(origin !== undefined)
        {
            var tr1 = new L.Matrix();
            tr1.m[3][1] = -origin.x;
            tr1.m[3][2] = -origin.y;
            rotMat.mul(tr1);

            var tr2 = new L.Matrix();
            tr2.m[1][1] = cosR;
            tr2.m[1][2] = sinR;
            tr2.m[2][1] = -sinR;
            tr2.m[2][2] = cosR;
            rotMat.mul(tr2);

            var tr3 = new L.Matrix();
            tr3.m[3][1] = origin.x;
            tr3.m[3][2] = origin.y;
            rotMat.mul(tr3);
        }
        else
        {
            rotMat.m[1][1] = cosR;
            rotMat.m[1][2] = sinR;
            rotMat.m[2][1] = -sinR;
            rotMat.m[2][2] = cosR;
        }
        this.m.mul(rotMat);
        return this;
    },
    
    getScale: function(){
        return [this.m.m[1][1], this.m.m[2][2]];
    },

    resetScale: function(){
        this.m.m[1][1] = 1;
        this.m.m[2][2] = 1;

        return this;
    },

    mapPoint: function(p) {
        var rx = p.x * this.m.m[1][1] + p.y * this.m.m[2][1] + this.m.m[3][1];
        var ry = p.x * this.m.m[1][2] + p.y * this.m.m[2][2] + this.m.m[3][2];
        p.x = rx;
        p.y = ry;
        return p;
    },

    toString: function(){
        return 'matrix(' +
                this.m.m[1][1] +', ' +// A
                this.m.m[2][1] +', ' +// C
                this.m.m[1][2] +', ' +// B
                this.m.m[2][2] +', ' +// D
                this.m.m[3][1] +', ' +// E = tX
                this.m.m[3][2] +      // F = tY
                ')';
    },

    toString3D: function(){
        /*
         * [A,B,0,0
         *  C,D,0,0
         *  E,F,1,0
         *  0,0,0,1]
         */
        return 'matrix3d(' +
                this.m.m[1][1] +','+ this.m.m[2][1] +',0,0,' +
                this.m.m[1][2] +','+ this.m.m[2][2] +',0,0,' +
                '0,0,1,0,' +
                this.m.m[3][1] +','+ this.m.m[3][2] +',0,1)';
    }
};


L.Transformation = function (a, b, c, d) {
	this._a = a;
	this._b = b;
	this._c = c;
	this._d = d;
};

L.Transformation.prototype = {
	transform: function (point, scale) { // (Point, Number) -> Point
		return this._transform(point.clone(), scale);
	},

	// destructive transform (faster)
	_transform: function (point, scale) {
		scale = scale || 1;
		point.x = scale * (this._a * point.x + this._b);
		point.y = scale * (this._c * point.y + this._d);
		return point;
	},

	untransform: function (point, scale) {
		scale = scale || 1;
		return new L.Point(
		        (point.x / scale - this._b) / this._a,
		        (point.y / scale - this._d) / this._c);
	}
};
