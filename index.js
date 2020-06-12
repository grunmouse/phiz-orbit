	var 
		math = Math,
		pi = math.PI,
		twoPi = 2*pi,
		sin = math.sin,
		cos = math.cos,
		tan = math.tan,
		acos = math.acos,
		asin = math.asin,
		sqrt = math.sqrt,
		sq = function(x){return x*x;},
		pow = math.pow,
		abs = math.abs,
		min = math.min,
		max = math.max,
		log = math.log,
		piphagorus = math.hypot,
		atan2 = math.atan2,
		angle = function(x, y){
			return atan2(y, x);
		},
		sgn=math.sign, //сигнум
		asinh=math.asinh, //гиперболический арксинус
		acosh=math.acosh, //гиперболический арккосинус
		round = math.round,
		
		Gamma = require('@grunmouse/phiz-const').GAMMA,
		
		Vector = require('@grunmouse/math-vector').Vector3,
		//Vector = linear.Vector,
		//Basis = linear.Basis,
		{ortRadius, ortTransversal, ortNormal, ortK} = require('./ort.js'),
		geometry = require('./geometry.js'),
		Radius = require('./radius.js'),
		{Period} = require('./period.js'),
		withS = require('./s.js');
		

	
	/**
	 * Интеграл энергии
	 * @param mu - гравитационный параметр
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 */
	function Energy(mu, p, eps){
		return -mu/p*(1-eps*eps);
	}
	
	/**
	 * Интеграл энергии через вектор состояния
	 * @param mu - гравитационный параметр
	 * @param r - радиус-вектор
	 * @param v - вектор скорости
	 */
	function EnergyV(mu, r, v){
		return v.eq - 2*mu/r.abs;
	}
	
	/**
	 * Перицентральное время для заданного радиуса
	 * симметрично относительно перицентра, неопределено для круговой орбиты
	 * @param mu - гравитационный параметр
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param r - радиус
	 */
	function t_r(mu, p, eps, r){
		var z, t, a, eps2, scale = 24*3600;
		
		mu = mu * scale*scale;
		
		if(eps == 1){
			t = (p+r)/sqrt(mu)*sqrt(2*r - p)/3;//Так сказал Wolfram
		}
		else{
			eps2 = eps*eps;
			a = p/abs(eps2-1) /* м */;
			z = parseFloat((1/eps+r*(eps2-1)/p/eps).toFixed(13));
			t = sqrt(a/mu)*a;
			if(eps<1){
				if(z<-1){
					z=-1;
				}
				else if(z>1){
					z=1;
				}
				t = t*(acos(z)-eps*sqrt(1-z*z));
			}
			else if(eps>1){
				if(z>-1 && z<1){
					z=1;
				}
				t = t*(eps*sqrt(z*z-1)-acosh(z));
			}
		}
		return t * scale;
	};
	
	/**
	 * Перицентральное время для заданного угла
	 * @param mu - гравитационный параметр
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function t_fi(mu, p, eps, fi){
		//fi = n*twoPi+fi_min
		//t = n*T + sgn(fi_min)*t(fi_min)
		var me = this,
			fi_min=fi, s, n=0, t=0, r;
		for(; fi_min>pi; fi_min-=twoPi){
			++n;
		}
		for(; fi_min<-pi; fi_min+=twoPi){
			--n;
		}
		if(eps<1){
			if(n){
				t = n*Period(mu, p, eps);
			}
			r = Radius(p, eps, fi_min);
			t += sgn(fi_min)*t_r(mu, p, eps, r);
		}
		else if(n==0){
			r = Radius(p, eps, fi_min);
			t = sgn(fi_min)*t_r(mu, p, eps, r);
		}
		else{
			return sgn(fi)*Infinity;
		}
		return t;
	}
	
	
	
	/**
	 * Радиальная составляющая скорости
	 * @param mu - гравитационный параметр
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function Vradial(mu, p, eps, fi){
		return sqrt(mu/p)*eps*sin(fi);
	}
	
	/**
	 * Трансверсальная составляющая скорости
	 * @param mu - гравитационный параметр
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function Vtransversal(mu, p, eps, fi){
		return sqrt(mu/p)*(1+eps*cos(fi));
	}
	
	/**
	 * Радиус-вектор тела на орбите
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param om - аргумент перицентра
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function vekRadius(Om, i, om, p, eps, fi){
		var u = om+fi;
		return ortRadius(Om, i, u).mul(Radius(p, eps, fi));
	}	
	
	/**
	 *	Орт вектора Лапласа равен орту радиус-вектора при fi=0
	 */
	function ortLaplas(Om, i, om){
		return ortRadius(Om, i, om)
	}
	
	/**
	 * Вектор скорости тела на орбите
	 * @param mu - гравитационный параметр
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param om - аргумент перицентра
	 * @param p - параметр орбиты
	 * @param eps - эксцентриситет орбиты
	 * @param fi - перицентральный угол (истинная аномалия)
	 */
	function vekV(mu, Om, i, om, p, eps, fi){
		var u = om+fi;
		return Vector.add(
			ortRadius(Om, i, u).mul(Vradial(mu, p, eps, fi)),
			ortTransversal(Om, i, u).mul(Vtransversal(mu, p, eps, fi))
		);
	}
	
	function vekLaplas(mu, Om, i, om, p, eps){
		var r = vekRadius(Om, i, om, p, eps, 0),
			v = vekV(mu, Om, i, om, p, eps, 0),
			c = r.cross(v); //интеграл площадей
		return r.ort().mul(-mu).add(v.cross(c)); // -mu*r.ort + (v x c) вектор Лапласа (направлен в перицентр)
	}
	
	/**
	 * Параметры орбиты по вектору состояния
	 * @param mu - гравитационный параметр
	 * @param r - радиус-вектор
	 * @param v - вектор скорости
	 */
	function mapOrbit(mu, r, v){

		var c = r.cross(v), //интеграл площадей
			i = acos(c.cosG), // =c.z/c.abs
			Om = atan2(c.x, -c.y),
			k = ortK(Om),
			l = r.ort().mul(-mu).add(v.cross(c)), // -mu*r.ort + (v x c) вектор Лапласа (направлен в перицентр)
			om = atan2(l.cross(k).abs, l.dot(k)),
			p = c.sq/mu,
			eps = l.abs/mu,
			fi = Vector.relAngle(l, r, c); //fi >0 , если (l, r, c) - правая тройка
			//cosfi = 1/eps*(p/r.abs - 1) - проверить
			
		/**
		 * @param Om - долгота восходящего узла
		 * @param i - наклонение
		 * @param om - аргумент перицентра
		 * @param p - параметр орбиты
		 * @param eps - эксцентриситет орбиты
		 * @param fi - перицентральный угол (истинная аномалия)
		 */
		return {
			Om:Om,
			i:i,
			om:om,
			p:p,
			eps:eps,
			fi:fi,
			c:c
		};
	}
	
	
	/**
	 * Проверяет вектор на компланарность с плоскостью орбиты
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param vek - проверяемый вектор
	 */
	function isInPlant(Om, i, vek){
		var k = ortK(Om),
			n = k.cross(vek);
		return n.abs==0 || abs(n.cosG) == abs(cos(i));
	}
	
	function basisOrbit(Om, i){
		var z = ortNormal(Om, i),
			x = ortK(Om),
			y = z.cross(x);
		return [x, y, z];
	}
	
	/**
	 *	Находит базис, в котором орбита описывается формулами x = R(fi)*cos(fi), y = R(fi)*sin(fi), где R - полярное уравнение конического сечения
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param om - аргумент перицентра
	 *	
	 */
	function basisCone(Om, i, om){
		var r = vekRadius(Om, i, om, 1, 0, 0),
			v = vekV(1, Om, i, om, 1, 0, 0),
			c = r.cross(v);	
		var z = c.ort(), //ortNormal(Om, i),
			x = r.ort(), //ortLaplas(Om, i, om),
			y = v.ort().neg(); //z.cross(x).neg();
			
		return [x, y, z];
	}
	
	/**
	 * Аргумент широты для радиус вектора
	 * @param Om - долгота восходящего узла
	 * @param i - наклонение
	 * @param r - радиус-вектор
	 */
	function calcUWithOrbit(Om, i, r){
		var k = ortK(Om),
			n = ortNormal(Om, i);
		return Vector.relAngle(k, r, n);
	}
	
	/**
	 *	Приблизительный радиус сферы действия планеты
	 *	@param mu0 - гравитационный параметр центрального тела
	 *	@param mu1 - гравитационный параметр планеты
	 *	@param ro - длина радиус-вектора планеты
	 */
	function Rs(mu0, mu1, ro){
		return ro*pow(mu1/mu0, 0.4);
	}
	
	
	/**
	 *	ХС перехода с круговой орбиты на гиперболическую
	 *	@parem v0 - скорость на круговой орбите
	 *	@parem vH - необходимый гиперболический остаток скорости
	 */
	function dVhyp(v0, vH){
		return v0*(sqrt(2+sq(vH/v0)) - 1);
	}
	
	module.exports = {
		/**
		 * Уравнение конического сечения в полярных координатах
		 * @param p - параметр орбиты
		 * @param eps - эксцентриситет орбиты
		 * @param fi - перицентральный угол (истинная аномалия)
		 */
		Radius:Radius,
		Period:Period,
		IntegralEnergy:Energy,
		IntegralEnergyV:EnergyV,
		vekRadius:vekRadius,
		vekV:vekV,
		V:{
			radial:Vradial,
			transrersal:Vtransversal
		},
		mapOrbit:mapOrbit,
		withS,
		tau:{
			r:t_r,
			fi:t_fi,
			s:withS.tau
		},
		isInPlant:isInPlant,
		calcUWithOrbit:calcUWithOrbit,
		ortK:ortK,
		ortNormal:ortNormal,
		ortLaplas:ortLaplas,
		basisCone:basisCone,
		Rs:Rs,
		geometry
		
	};