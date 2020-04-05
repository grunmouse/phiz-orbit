let {PI, sqrt, asin, abs, cos, sin, atan2} = Math;
let common = require('./common.js');
let {sin2a, cos2a} = require('./trigonometr.js');
let {ellipse} = require('./f_g_equation.js');
let NewRaf = require('./newton_raffson.js');


function iterFi(mu, r1, r2, tau, s, eps){
	let {dfi_cos, dfi_sin} = common(r1, r2, s);
	dfi_sin = dfi_sin*s;
	
	let dfi = atan2(dfi_sin, dfi_cos);
	
	let fi1 = 0;
	while(calcAE(fi1, dfi, r1.abs, r2.abs).e<0){
		fi1+= PI/18; //10°
	}
	while(calcAE(fi1, dfi, r1.abs, r2.abs).a<0){
		fi1+= PI/18; //10°
	}
	
	let F = (fi1)=>iter(fi1, dfi, mu, r1.abs, r2.abs, tau);
	let dF = (fi1)=>((F(fi1+PI/100) - F(fi1))/(PI/100));
	fi1 = NewRaf(F, dF, fi1, eps);
	
	let {fi2, e, a, E1, E2} = calc(fi1, dfi, mu, r1.abs, r2.abs);
	
	let f = ellipse.f(a, r1, E2-E1);
	
	let g = ellipse.g(mu, a, E2-E1, tau);
	
	let v1 = r2.add(r1.mul(-f)).mul(1/g);

	return v1;
	
}

function calcAE(fi1, dfi, r1, r2){
	let fi2 = fi1+dfi;
	let e = (r2-r1)/(r1*cos(fi1)-r2*cos(fi2));
	let a;
	if(e>=0){
		a = r1*(1+e*cos(fi1))/(1-e**2);
	}
	return {fi2, e, a};
}

function calc(fi1, dfi, mu, r1, r2){
	let {fi2, e, a} = calcAE(fi1, dfi, r1, r2);
	let E1_sin, E1_cos, E2_sin, E2_cos, E1, E2, n;
	if(e>=0){
		a = r1*(1+e*cos(fi1))/(1-e**2);
		if(a>0){
			E1_sin = sqrt(1-e**2)*sin(fi1)/(1+e*cos(fi1));
			E1_cos = (cos(fi1) + e)/(1+e*cos(fi1));
			E2_sin = sqrt(1-e**2)*sin(fi2)/(1+e*cos(fi2));
			E2_cos = (cos(fi2) + e)/(1+e*cos(fi2));
			E1 = atan2(E1_sin, E1_cos);
			E2 = atan2(E2_sin, E2_cos);
			n = sqrt(mu/a**3);
		}
	}
	return {e, a, E1, E2, E1_sin, E2_sin, n};
}

function iter(fi1, dfi, mu, r1, r2, tau){
	let {e, a, E1, E2, E1_sin, E2_sin, n} = calc(fi1, dfi, mu, r1, r2);

	let M1 = E1 - e*E1_sin;
	let M2 = E2 - e*E2_sin;
	let F = tau - (M2-M1)/n;
	return F;
	
}

module.exports = iterFi;