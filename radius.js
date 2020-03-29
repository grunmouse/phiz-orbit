	/**
	 * Уравнение конического сечения в полярных координатах
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function Radius(p, eps, fi){
		return p/(1+eps*Math.cos(fi));
	}
	
module.exports = Radius;