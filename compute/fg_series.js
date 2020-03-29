let recursiveFunction = require('@rakov/recursive-function.js');

let F = recursiveFunction([1], (F, n)=>(F(n-1)*n));

function C(n, k){
	return F(n)/F(k)/F(n-k);
}

/**
 * Суммирует массивы а и b перемножая их попарно и домножая на биномиальный коэффициент
 * @param a - массив
 * @param b - массив
 * @param n - степень
 */
function binom(a, b, n){
	let sum = 0;
	for(let k=0; k<=n; ++k){
		sum += C(n,k)*a[n-k]*b[k];
	}
	return sum;
}

function tailor(f, x, n){
	let sum = 0;
	for(let k=0; k<=n; ++k){
		sum += f[k]*x**k /F(k);
	}
	return sum;
}

function getScalar(mu, vecR, dvecR){
	let r = vecR.abs, V = dvecR.abs, V_sq = dvecR.sq, r_sq = vecR.sq;
	let u = mu/r**3;
	let p = vecR.smul(dvecR)/r_sq;
	let q = V_sq/r_sq - u;
	return {r, V, u, p, q};
	
}

function getScalar2(mu, r0, dr0, V0){
	let r = r0, V = V0, V_sq = V**2, r_sq = r**2;
	let u = mu/r**3;
	let p = dr0/r0;
	let q = V_sq/r_sq - u;
	return {r, V, u, p, q};
	
}

function getUPQ(mu, r0, dr0, V0, n){
	let sc = getScalar2(mu, r0, dr0, V0);
	let u=[sc.u], p=[sc.p], q=[sc.q];
	u[1] = -3*u[0]*p[0];
	p[1] = q[0] - 2*p[0]**2;
	q[1] = -p[0]*(u[0] + 2*q[0]);
	for(let k=2; k<=n; ++k){
		u[k] = -3*binom(u, p, k-1);
		p[k] = q[k-1] - 2*binom(p, p, k-1);
		q[k] = -binom(p, u.map((ui, i)=>(ui + 2*q[i])), k-1);
	}
	return {u, p, q};
}

function f0(u, n){
	let f = [1, 0, -u[0], -u[1]];
	for(let i = 4; i<=n; ++i){
		let m = i-2;
		let sum = u[m];
		for(let k=2; k<=m; ++k){
			sum += C(n-2, k)*u[n-2-k]*f[k];
		}
		f[i] = -sum;
	}
	return f;
}

function g0(u, n){
	let g = [0, 1, 0, -u[0]];
	for(let i = 4; i<=n; ++i){
		let m = i-2;
		let sum = C(n-2, 1)* u[i-3];
		for(let k=2; k<=m; ++k){
			sum += C(n-2, k)*u[n-2-k]*g[k];
		}
		g[i] = -sum;
	}
	return g;
}

function fg(mu, r0, dr0, V0, tau, n){
	let {u} = getUPQ(mu, r0, dr0, V0, n-2);
	let f = f0(u, n), g = g0(u, n);
	
	return {
		f:tailor(f, tau, n);,
		g:tailor(g, tau, n)
	};
}

function dfg(mu, r0, dr0, V0, tau, n){
	let {u} = getUPQ(mu, r0, dr0, V0, n-2);
	let f = f0(u, n), g = g0(u, n);
	f.shift(); g.shift();
	return {
		df:tailor(f, tau, n);,
		dg:tailor(g, tau, n)
	};
}

function f(V0, r0, dr0, tau){
}

function g(V0, r0, dr0, tau){

}

function df(V0, r0, dr0, tau){

}

function dg(V0, r0, dr0, tau){

}

module.exports = {
	fg, dfg
}