let Vector = require('@grunmouse/math-vector').Vector3;
let common = require('./common.js');
let {sqrt, acos, sign, abs} = Math;
let {ellipse} = require('./f_g_equation.js');

function gauss(mu, r1, r2, tau, s, eps){
	let {dfi_cos} = common(r1, r2, s);
	let dfi_2_div_cos = sqrt((1+dfi_cos)/2);
	
	let l = (r1.abs+r2.abs)/(4 *sqrt(r1.abs*r2.abs)*dfi_2_div_cos) -0.5;
	
	let m = mu*tau**2/(2*sqrt(r1.abs*r2.abs)*dfi_2_div_cos)**3;
	
	let {y, dE_sub_2_div_sin, dE_sub, dE_sub_sin, dE_sub_cos} = iter(m, l, eps);
	
	let a_sqrt = tau*sqrt(mu)/(2*y*sqrt(r1.abs*r2.abs)*dfi_2_div_cos*dE_sub_2_div_sin);
	
	let a = a_sqrt**2;
	
	let f = ellipse.f(a, r1.abs, dE_sub);
	
	let g = ellipse.g(mu, a, dE_sub, tau);
	
	let v1 = r2.add(r1.mul(-f)).mul(1/g);
	
	return v1;
}

function iter(m, l, eps){
	
	let y=1; //Первое приближение
	let dE_sub_2_div_sin, dE_sub, dE_sub_sin, dE_sub_cos, y_old;
	do{
		let x = m/y**2 - l;
		let dE_sub_2_div_cos = 1-2*x;
		
		dE_sub_2_div_sin = sqrt(4*x*(1-x));
		
		dE_sub_sin = 2 * dE_sub_2_div_sin * dE_sub_2_div_cos;
		
		dE_sub_cos = dE_sub_2_div_cos**2 - dE_sub_2_div_sin**2;
		
		dE_sub = acos(dE_sub_cos)*sign(dE_sub_sin);
		
		let X = (dE_sub - dE_sub_sin)/dE_sub_2_div_sin**3;
		
		y_old = y;
		y = 1+X*(l+x);
	}while(abs(y_old-y)>eps);
	return {y, dE_sub_2_div_sin, dE_sub, dE_sub_sin, dE_sub_cos};
}

module.exports = gauss;