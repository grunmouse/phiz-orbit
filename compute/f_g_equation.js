let {sqrt, sin, cos, cosh, sinh} = Math;
/**
 * @param a - большая полуось
 * @param r1 - модуль первого радиус-вектора
 * @param dEE - E2-E1 - разность эксцентрических аномалий второго и первого радиус-векторов
 */
function ellipse_f(a, r1, dEE){
	return 1 - a/r1*(1 - cos(dEE));
}
/**
 * @param mu - гравитационный параметр центрального тела
 * @param a - большая полуось
 * @param dEE - E2-E1 - разность эксцентрических аномалий второго и первого радиус-векторов
 * @param tau - время полёта из точки 1 в току 2
 */
function ellipse_g(mu, a, dEE, tau){
	return tau - sqrt(a**3/mu)*(dEE - sin(dEE));
}

function ellipse_df(mu, a, r0, r, dEE){
	return -sqrt(mu*a)/r/r0*sin(dEE);
}

function ellipse_dg(a, r, dEE, Ce, Se){
	return a/r*(cos(dEE) - Ce*cos(dEE) + Se*sin(dEE));
}

function hyperbola_f(a, r0, dFF){
	return 1 - a/r0*(1-cosh(dFF));
}

function hyperbola_g(mu, a, dFF, tau){
	return tau - sqrt(-(a**3)/mu)*(sinh(dFF) - dFF);
}

module.exports = {
	ellipse:{
		f:ellipse_f,
		g:ellipse_g,
		df:ellipse_df,
		dg:ellipse_dg
	},
	hyperbola:{
		f:hyperbola_f,
		g:hyperbola_g
	}
}