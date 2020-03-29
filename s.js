	var 
		Stumpf = require('@grunmouse/math-stumpff'),
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
		round = math.round;

	/**
	 * время полёта из точки, заданной радиус-вектором, до точки, заданной фиктивным временем
	 * @param mu - гравитационный параметр
	 * @param r - радиус-вектор
	 * @param v - вектор скорости
	 * @param s - фиктивное время
	 */
	function tau(mu, r, v, s){
		var C = Stumpf(xStumpf(mu, r, v, s));
		
		return s*(r.abs*C[1] + s*(r.smul(v)*C[2] + s*mu*C[3]));	//r*s*C1 + (r·v)*(s**2)*C2 + mu*(s**3)*C3
	}
	
	/**
	 *	Аргумент x для функций Штумпфа
	 */
	function xStumpf(mu, r, v, s){
		var h = EnergyV(mu, r, v),
			x = -h*s*s;
			
		return x;
	}
	
	/**
	 *	Расстояние от центрального тела через фиктивное время s
	 */
	function r_s(mu, r, v, s){
		var C = Stumpf(xStumpf(mu, r, v, s));
		
		return r.abs*C[0] + s*(r.smul(v)*C[1] + mu*s*C[2]); // r*C0 + (r·v)*s*C1 + mu*(s**2)*C2
	}
	
	/**
	 *	вычитаемое следующего приближения фиктивного времени для заданного tau
	 */
	function ds_next(mu, absR, rv, h, tau, s){
		var C = Stumpf(-h*s*s);
		
		var smu = s*mu;
		
		return (s*(absR*C[1] + s*(rv*C[2] + smu*C[3]))-tau)/(absR*C[0] + s*(rv*C[1] + smu*C[2]));
	}
	
	/**
	 *	Начальное приближение фиктивного времени для заданного tau
	 */
	function s_0(mu, h, tau){
		return h<0 ? -h*tau/mu : 0;
	}
	
	/**
	 *	Численное вычисление фиктивного времени для заданных начального вектора состояния и времени полёта
	 *	@param mu - гравитационный параметр
	 *	@param r - начальный радиус-вектор
	 *	@param v - начальная скорость
	 *	@param tau - время полёта
	 *	@param eps - допустимая погрешность результата
	 */
	function calcS(mu, r, v, tau, eps){
		var h = EnergyV(mu, r, v);	
		var rv = r.smul(v);
		var absR = r.abs;
		
		var s = s_0(mu, h, tau);
		var ds = ds_next(mu, absR, rv, h, tau, s);
		while(ds>eps){
			s = s - ds;
			ds = ds_next(mu, absR, rv, h, tau, s);
		}
		return s;
	}
	
	function vekRadius_s(mu, r, v, tau, s){
		var C = Stumpf(xStumpf(mu, r, v, s));
		var mus2 = mu*s*s;
		var f = 1 - mus2*C[2]/r,
			g = tau - mus2*s*C[3];
			
		return r.mul(f).sum(v.mul(g));
	}
	
	function vekV_s(mu, r, v, s){
		var C = Stumpf(xStumpf(mu, r, v, s));
		var mus = mu*s;
		var R = r.abs*C[0] + s*(r.smul(v)*C[1] + mus*C[2]);
		var mus_r = mus/R;
		var 
			df = - mus_r*C[1]/R,	//df/ds = - mu*s*C1/(r*R)
			dg = 1 - mus_r*s*C[2];	//dg/ds = 1 - mu*(s**2)*C2/R
		return r.mul(df).sum(v.mul(dg));
	}
	
	function Lambert(mu, r0, r1, tau, eps){
		let rsum = (r0.abs + r1.abs);
		
		let rsum_mu_div_sqrt = sqrt(rsum/mu)
		
		let fi = r0.angle(r1);
		let ro = sqrt(2*r0.abs*r1.abs)/rsum*cos(fi/2);
		let sigma = sqrt(mu/rsum**3)*tau;
		let sigma_par = 1/3*(Math.SQRT2+ro)*sqrt(1-sqrt(2*ro));
		let x_min = -4*log(1/Math.SQRT2/ro * sqrt(0.5/ro**2 - 1))**2;
		
		function F(x){
			let C = Stumpf(x, 7);
			let C2_sqrt = sqrt(C[2]);
			let C2_3_pow_sqrt = sqrt(C[2]**3);
			
			let u = sqrt(1-ro*C[1]/C2_sqrt);
			
			let u_C2_sqrt_div = u/C2_sqrt;
			
			let u_C2_sqrt_div_pow_3 = u_C2_sqrt_div**3;
			let s = rsum_mu_div_sqrt * u_C2_sqrt_div;
			
			let u_3_pow = u**3;
			let value = C[3]/C2_3_pow_sqrt*u_3_pow + ro*u;
			
			let dF_dx = (C[3]**2 - C[5] + 4*C[6])/4/C2_3_pow_sqrt*u_3_pow +(3*C[3]*u**2/C2_3_pow_sqrt + ro)*ro*C2_sqrt/8/u;
			//sukhanov 2010 str 69
			return {
				F:value,
				dF:dF_dx
			}
		}
		function dx(x){
			let {F, dF} = F(x);
			return (F-sigma)/Fx;
		}
		if(sigma<=sigma_par){
			//hiperbola
			//parabola
			x=Math.max(0, x_min);
		}
		else if(sigma>sigma_par){
			//ellipse
			let um = sqrt(1-Math.SQRT2*abs(ro));
			let e0 = (PI/(2*um**3/3+sigma-ro*um))**(1/3)*um;
			x = 4*(PI-e0)**2;
		}
		
		do{
			let d = dx(x);
			x = x - d;
		}while(d>eps);
		
		
	}

module.exports = {
	tau,
	v:vekV_s,
	r:vekRadius_s,
	calcS,
	Lambert
}