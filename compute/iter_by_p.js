let {sqrt, asin, abs, cos, sin, atan2} = Math;
let common = require('./common.js');
let {sin2a, cos2a} = require('@rakov/math-trigonometr');
let {ellipse} = require('./f_g_equation.js');
let NewRaf = require('./newton_raffson.js');

function iterP(mu, r1, r2, tau, s, eps){
	let {dfi_cos, dfi_sin} = common(r1, r2, s);
	dfi_sin = dfi_sin*s;
	let p = (r1.abs+r2.abs)*0.4;
	
	let F = (p)=>iter(p, r1.abs, r2.abs, dfi_cos, dfi_sin);
	let dF = (p)=>((F(1.05*p)-F(p))/(p*0.05);
	
	p = NewRaf(F, dF, p, eps);
	
	let {a, E1, E2, E1_sin, E2_sin} = calc(p, mu, r1.abs, r2.abs, dfi_cos, dfi_sin);
	
	let f = ellipse.f(a, r1, E2-E1);
	
	let g = ellipse.g(mu, a, E2-E1, tau);
	
	let v1 = r2.add(r1.mul(-f)).mul(1/g);

	
	return v1;
	
}

function calc(p, mu, r1, r2, dfi_cos, dfi_sin){
	let e_fi1_cos_mul = p/r1 - 1;
	let e_fi2_cos_mul = p/r2 - 1;
	let e_fi1_sin_mul = (e_fi1_cos_mul*dficos - e_fi2_cos_mul)/dfi_sin;
	let e_fi2_sin_mul = (-e_fi2_cos_mul*dficos + e_fi1_cos_mul)/dfi_sin;
	let e_sq = e_fi1_sin_mul**2 + e_fi1_cos_mul**2;
	let a = p/(1-e_sq);
	let n = sqrt(mu/a**3);
	let E1_cos, E2_cos, E1_sin, E2_sin;
	let fi1_sin, fi1_cos, fi2_sin, fi2_cos, e;
	if(e!=0){
		e = sqrt(e_sq);
		fi1_sin = e_fi1_sin_mul/e;
		fi2_sin = e_fi2_sin_mul/e;
		fi1_cos = e_fi1_cos_mul/e;
		fi2_cos = e_fi2_cos_mul/e;
		E1_cos = r1/p*(fi1_cos+e);
		E2_cos = r2/p*(fi2_cos+e);
		E1_sin = r1/p*sqrt(1-e_sq)*fi1_sin;
		E2_sin = r2/p*sqrt(1-e_sq)*fi2_sin;
	}
	else{
		e =0;
		E1_cos = 1;
		E2_cos = dfi_cos;
		E1_sin = 0;
		E2_sin = dfi_sin;
	}
	let E1 = atan2(E1_sin, E1_cos);
	let E2 = atan2(E2_sin, E2_cos);
	return {a, E1, E2, E1_sin, E2_sin};

}

function iter(p, r1, r2, dfi_cos, dfi_sin){
	let {E1, E2, E1_sin, E2_sin} = calc(p, r1, r2, dfi_cos, dfi_sin);
	
	let M1 = E1 - e*E1_sin;
	let M2 = E2 - e*E2_sin;
	let F = tau - (M2-M1)/n;
	return F;
}