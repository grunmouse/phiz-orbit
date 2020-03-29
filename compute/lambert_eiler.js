let {sqrt, asin, abs, cos, sin} = Math;
let common = require('./common.js');
let {sin2a, cos2a} = require('@rakov/math-trigonometr');
let {ellipse} = require('./f_g_equation.js');
let NewRaf = require('./newton_raffson.js');

function lambert_eiler(mu, r1, r2, tau, s, eps){
	let {dfi_cos, dfi_sin} = common(r1, r2, s);
	let dfi_2_div_cos = sqrt((1+dfi_cos)/2);
	dfi_sin *= s;
	let a = (r1.abs + r2.abs)/2;
	let c = sqrt(r1.sq + r2.sq - 2 *r1.smul(r2));
	
	let Fa  = (a)=>funcF(mu, r1.abs, r2.abs, tau, dfi_2_div_cos, a, c, s);
	let da = (a)=>(a*0.05);
	let dFa = (a)=>((Fa(a*1.05)-Fa(a))/(a*0.05));
	
	a = NewRaf(Fa, dFa, a, eps);
	
	let {e, del, e_sin, del_sin} = e_del(r1.abs, r2.abs, dfi_2_div_cos, a, c, s)
	
	let E2_E1_sub = e - del;
	
	let f = ellipse.f(a, r1, E2_E1_sub);
	
	let g = ellipse.g(mu, a, E2_E1_sub, tau);
	
	let v1 = r2.add(r1.mul(-f)).mul(1/g);

	
	return v1;
}

function e_del(r1, r2, dfi_2_div_cos, a, c, s){
	let e_2_div_sin = sqrt((r2 + r1 +c)/4/a);
	let del_2_div_sin = sqrt(r2*r1)*dfi_2_div_cos/2/a/e_2_div_sin;
	let del_2_div_cos = sqrt(1-(r2 + r1 - c)/4/a);
	let e_2_div_cos = s*sqrt(1-e_2_div_sin**2);
	let del_sin = sin2a(del_2_div_sin, del_2_div_cos);
	let e_sin = sin2a(e_2_div_sin, e_2_div_cos);
	let e = asin(e_sin), del = asin(del_sin);
	return {e, del, e_sin, del_sin};
}

function funcF(mu, r1, r2, tau, dfi_2_div_cos, a, c, s){
	let {e, del, e_sin, del_sin} = e_del(r1, r2, dfi_2_div_cos, a, c, s);
	
	return tau - sqrt(a**3/mu)*((e-e_sin) - (del -del_sin));
	
}

module.exports = lambert_eiler;