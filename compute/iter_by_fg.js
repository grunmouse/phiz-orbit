let {PI, sqrt, asin, abs, cos, sin, atan2} = Math;
let common = require('./common.js');
let {sin2a, cos2a} = require('@rakov/math-trigonometr');
let {ellipse} = require('./f_g_equation.js');
let NewRaf = require('./newton_raffson.js');

let series = require('./fg_series.js');

function iterFG(mu, r1, r2, tau, s, eps){
	let {dfi_cos, dfi_sin} = common(r1, r2, s);
	dfi_sin = dfi_sin*s;
	
	let dfi = atan2(dfi_sin, dfi_cos);
	
	//v_2 - v_1 < 90° , иначе процесс не сходящийся
	
	let t0 = tau/2;
	let tau1 = - t0;
	let tau2 = tau - t0;
	
	let r0 = (r2.abs+r1.abs)/2;
	let A = 1- mu*tau1**2/2/r0**3;
	let B = 1- mu*tau2**2/2/r0**3;
	let Delta = A*tau2 - B*tau1;
	let vecR0 = r1.mul(tau2/Delta).add(r2.mul(-tau1/Delta));
	let dvecR0 = r2.mul(A/Delta).add(r1.mul(-B/Delta));
	r0 = vecR0.abs;
	let V0, dr0;
	{r0, V0, dr0} = calcR(vecR0, dvecR0);
	
	let del_r0, del_V0, del_dr0;
	do{
		{vecR0, dvecR0} = iter(V0, r0, dr0, tau1, tau2, vecR1, vecR2);
		let {nr0, nV0, ndr0} = calcR(vecR0, dvecR0);
		del_V0 = abs(nV0-V0);
		del_dr0 = abs(ndr0-dr0);
		del_r0 = abs(nr0-r0);
		{r0, V0, dr0} = {nr0, nV0, ndr0};
	}while(del_r0>eps || del_V0>eps || del_dr0>eps); //TODO надо разные eps
	
	let d = dfg(mu, r0, dr0, V0, tau1, n)
	
	let v1 = vecR0.mul(d.df).add(dvecR0.mul(d.dg));
	
	return v1;
}

function calcR(vecR0, dvecR0){
	let r0 = vecR0.abs;
	let V0 = dvecR0.abs;
	let dr0 = vecR0.smul(dvecR0).mul(1/r0);
	return {r0, V0, dr0};
}

function calcC(V0, r0, dr0, tau1, tau2){
	let fg1 = series.fg(mu, r0, dr0, V0, tau1, 8);
	let fg2 = series.fg(mu, r0, dr0, V0, tau2, 8);
	let f1 = fg1.f, g1 = fg1.g;
	let f2 = fg2.f, g2 = fg2.g;
	let D = f1*g2 - f2*g1;
	let C1 = g2/D;
	let C2 = -g1/D;
	let dC1 = - f2/D;
	let dC2 = f1/D;
	
	return {C1, C2, dC1, dC2};
}

function iter(V0, r0, dr0, tau1, tau2, vecR1, vecR2){
	let {C1, C2, dC1, dC2} = calcC(V0, r0, dr0, tau1, tau2);
	let vecR0 = r1.mul(C1).add(r2.mul(C2));
	let dvecR0 = r1.mul(dC1).add(r2.mul(dC2));
	
	return {vecR0, dvecR0};	
}