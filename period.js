const {PI, sqrt} = Math;
const twoPI = 2*PI;

	/**
	 * Период обращения по эллиптической орбите
	 * @param mu - гравитационный параметр
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function Period(mu, p, eps){
		if(eps>=1) return Infinity;
		var a = p/(1-eps*eps);
		return twoPi*sqrt(a/mu)*a;
	}		
	
	
	function ellipseA(mu, T){
		/*
			T = twoPi*sqrt(a/mu)*a;
			sqrt(mu)*T/twoPi = a**(3/2);
			a = (sqrt(mu)*T/twoPi)**(2/3);
		*/
		
		return 
	}
		
	
module.exports = {
	Period
};