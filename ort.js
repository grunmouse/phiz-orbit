let Vector = require('@grunmouse/math-vector').Vector3;
let {cos, sin} = Math;

	/**
	 * Орт радиус-вектора тела на орбите
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param u - аргумент широты
	 */
	function ortRadius(Om, i, u){
		var cOm = cos(Om), 
			sOm = sin(Om), 
			ci = cos(i), 
			si = sin(i), 
			cu = cos(u),
			su = sin(u),
			suci = ci*su;
		return new Vector(
			cOm*cu - sOm*suci,
			sOm*cu + cOm*suci,
			su*si
		);
	}
	
	/**
	 * Орт трансверсального направления для тела на орбите
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param u - аргумент широты
	 */
	function ortTransversal(Om, i, u){
		var cOm = cos(Om), 
			sOm = sin(Om), 
			ci = cos(i), 
			si = sin(i), 
			cu = cos(u),
			su = sin(u),
			cuci = ci*cu;
		return new Vector(
			-cOm*su - sOm*cuci,
			-sOm*su + cOm*cuci,
			cu*si
		);
	}
	
	/**
	 * Орт нормали к орбите
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 */
	function ortNormal(Om, i){
		return new Vector(
			sin(Om)*sin(i),
			-cos(Om)*sin(i),
			cos(i)
		);
	}
	
	/**
	 * Орт направления на восходящий узел
	 * @param Om - долгота восходящего узла
	 */
	function ortK(Om){
		return new Vector(cos(Om), sin(Om), 0);
	}

module.exports = {ortRadius, ortTransversal, ortNormal, ortK};