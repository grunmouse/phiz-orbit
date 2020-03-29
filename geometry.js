var Radius = require('./radius.js');

var {genericA} = require('@grunmouse/generic-functions');

var sqrt = Math.sqrt;
//R = p/(1+eps*cos(fi))
//R(0) = p/(1+eps)
//R(PI) = p/(1-eps)

function r0(p, eps){
	return p/(1+eps);
}

function rPI(p, eps){
	if(eps<0){
		return p/(1-eps);
	}
	else{
		return Infinity;
	}
}

function p(r0, eps){
	return r0*(1+eps);
}

//Большая полуось
function ellipsA(p, eps){
	return p/(1-eps**2);
}

//Малая полуось
function ellipsB(p, eps){
	return p/sqrt(1-eps**2);
}

//Полуфокусное расстояние
function ellipsC(p, eps){
	return ellipsA(p, eps)*eps;
}

function ellipsFull(p, eps){
	let e = (1-eps**2);
	let a = p/e;
	return {
		a,
		b:p/sqrt(e),
		c:a*eps
	};
}

//Расстояние о центра до вершины
function hyperbolaA(p, eps){
	return p/(eps**2 - 1);
}

//Расстояние от фокуса до асимптоты (прицельный параметр)
function hyperbolaB(p, eps){
	return p/sqrt(eps**2-1);
}

//Расстояние от центра до фокуса
function hyperbolaC(p, eps){
	return hyperbolaA(p, eps)*eps;
}

function hyperbolaFull(p, eps){
	let e = (eps**2-1);
	let a = p/e;
	return {
		a,
		b:p/sqrt(e),
		c:a*eps
	};
}

//Для параболы - коэффициенты квадратичной функции x = a*y**2 + c;
function parabolaA(p){
	return -0.5/p;
}
function parabolaC(p){
	return p/2;
}
function parabolaFull(p){
	return {
		a:-0.5/p,
		c:p/2
	};
}

function fi(p, eps, r){
	return eps>0 ? Math.acos((p/r - 1)/eps) : NaN;
}

module.exports = {
	ellips:{
		A:ellipsA,
		B:ellipsB,
		C:ellipsC,
		full:ellipsFull
	},
	hyperbola:{
		A:hyperbolaA,
		B:hyperbolaB,
		C:hyperbolaC,
		full:hyperbolaFull
	},
	parabola:{
		A:parabolaA,
		C:parabolaC,
		full:parabolaFull
	},
	r0,
	rPI,
	fi,
	genP:genericA({
		"eps,r0":map=>p(map.r0, map.eps)
	},
	['r0', 'eps'])
}