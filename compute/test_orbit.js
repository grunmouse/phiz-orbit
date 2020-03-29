let orbit = require('../index.js');

//function form(


function twoR(mu, map, fi1, fi2){
	let {Om, i, om, p, eps} = map;
	return {
		tau:orbit.tau.fi(mu, p, eps, fi2) - orbit.tau.fi(mu, p, eps, fi1),
		r1:orbit.vekRadius(Om, i, om, p, eps, fi1),
		r2:orbit.vekRadius(Om, i, om, p, eps, fi2)
	}
}

let PI = Math.PI, twoPI = PI*2, {acos, sqrt, min} = Math;

/**
 * Генерирует пары {наклонение, долгота восходящего узла}
 * @rapam s - угловой шаг генерации
 */
function* rotateOrbits(s){
	yield {i:0, Om:0};
	for(let i = s; i<PI; i+=s){
		for(let Om = 0; Om<twoPI; om+=s){
			yield {i, Om};
		}
	}
	return;
}
/**
 * Генерирует пары {эксцентриситет, аргумент перицентра}
 * для эксцентриситетов 0...0.9 - если e===0
 *	или 1.1...1.9 - если e === 1;
 *	шаг эксцентриситета 0.1
 * @rapam s - угловой шаг генерации
 * @param e - начальный эксцентриситет 0 или 1, 1 не будет включена в результат. 
 */
function *epsSpinOrbit(s, e){
	if(e===0){
		yield {eps:0, om:0};
	}
	for(let eps = e+0.1; eps<e+1; e+=0.1){
		for(let om = 0; om<twoPI; om+=s){
			yield {eps, om};
		}
	}
	return;
}

/**
 * Генерирует орбиты для параметров, генерируемых itrP, c угловым шагом s и диапазоном эксцентриситетов e..e+0.9 не включая 1
 */
function *orbits(itrP, e, s){
	for(let {eps, om} of epsSpinOrbit(s, e)){
		for(let p of itrP){
			for(let {i, Om} of rotateOrbits(s)){
				yield {p, eps, i, Om, om};
			}
		}
	}
	return;
}

/**
 * Генерирует задачи распознавания орбит
 * @param mu - гравитационный параметр центрального тела
 * @param orbits - итерируемый набор орбит
 * @param f - угловой шаг набора радиус-векторов
 * @param df - наибольший угол между радиус-векторами
 * Генерируется пара векторов, такая, чтобы угол между ними был больше нуля, но меньше df
 *	для эллипса начальный угол не меньше -Pi, а конечный не больше 3*Pi
 *	для гиперболы - углы в пределах асимптот
 */
function *tasks(mu, orbits, f, df){
	df = df || twoPI;
	for(let map of orbits){
		let eps = map.eps;
		let teta = eps<1 ? PI : (acos(-sqrt((1-eps**2)/(2-eps**2)))-f);
		for(let fi1 = -teta; fi1<=teta-f; fi1+=f){
			let max2 = eps<1 ? fi1+df : min(fi1+df, teta);
			for(let fi2 = fi1+f; fi2<=max2; fi2+=f){
				let result = {...map, ...twoR(mu, map, fi1, fi2)};
				console.log(result);
				yield result;
			}
		}
	}
	console.log('ok');
	return;
}

module.exports = {
	twoR,
	orbits,
	tasks
}