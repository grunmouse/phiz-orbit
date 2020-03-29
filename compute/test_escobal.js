let {Vector} = require('@rakov/math-linear-d3');
/*

	k = 0.07436574 (э.р.)**(3/2)/мин
	\mu = 1 з.м. 
	a = 1 э.р. \approx 6371 км
	
	k**2 = 0.0055302632857476 э.р.**3/мин**2
	
	э.р. = 6371 км = 6371000 м
	
	\mu_s = k**2 * (м/э.р.)**3 * (мин/c)**2
	
	\mu_s = 3.9725202731797994e+14 м^3/c^2 - принято у Эскобала
	\mu_s = 3.9857605760000006e+14 м^3/c^2 - рассчитано по данным гугла
*/
let er = 6371000;
let mu =3.9725202731797994e+14; //как у Эскобала

function JD2sec(JD){
	return (JD - 2440587.5)*86400;
}
function dJD2sec(JD){
	return (JD)*86400;
}

function date2JD(...a){
	(Date.UTC(...a) / 86400000) + 2440587.5
}

function date2sec(...arg){
	return Date.UTC(...arg).valueOf()/1000;
}

let deg = Math.PI/180;



let sample1 = {
	tau: dJD2sec(0.01044412), //JD
	r1:new Vector(2.460809*er, 2.040523*er, 0.143819*er),
	r2:new Vector(1.988041*er, 2.503334*er, 0.314554*er),
	a:4*er,
	eps:0.2,
	i:15*deg,
	Om:30*deg,
	om:10*deg,
	T:date2sec(1964, 0, 01) // "1964-01-01 0:0:0"
};

let sample2 = {
	tau: dJD2sec(0.01527809),
	r1: new Vector(-1.75981, 1.68113, 1.16913).mul(er),
	r2: new Vector(-2.23077, 0.77454, 1.34602).mul(er),
	a:3*er,
	eps:0.1,
	i:30*deg,
	Om:80*deg,
	om:60*deg,
	T:date2sec(1964, 5, 01) // "1964-01-01 0:0:0"
};

let sample3 = {
	tau: dJD2sec(0.01316924),
	r1:new Vector(0.41136, -1.66250, 0.82272).mul(er),
	r2:new Vector(0.97757, -1.64428 /* ? */, -0.042363).mul(er),
	a:2*er,
	eps:0.05,
	i:60*deg,
	Om:120*deg,
	om:150*deg,
	T:date2sec(1963,11,23 )
}

let sample4 = {
	tsu:dJD2sec(0.14971172),
	r1:new Vector(2.78418, 0.82815, 0.75).mul(er),
	r2:new Vector(2.3728, 1.54778, 1.11792).mul(er),
	a:6*er,
	eps:0.5,
	i:150*deg,
	Om:170*deg,
	om:150*deg,
	T:date2sec(1964, 0, 15, 22, 30, 3)
}

let sample5 = {
	tau:dJD2sec(0.42710476),
	r1:new Vector(-1.00746, -3.89711, -2.01185).mul(er),
	r2:new Vector(-0.64633, 5.12188, -1.29069).mul(er),
	a:0.5*er,
	eps:0.1,
	i:63.4*deg,
	Om:270*deg,
	om:330*deg,
	T:date2sec(1964, 7,4, 12)
}

let sample6 ={
	tau:dJD2sec(0.21227310),
	r1:new Vector(-2.57823, 2.13649, 0.59004).mul(er),
	r2:new Vector(3.49838, -2.9461, 0.23276).mul(er),
	a:4*er,
	eps:0.15,
	i:88*deg,
	Om:140*deg,
	om:10*deg,
	T:date2sec(1964, 0, 30)
}

module.exports = {
	sample1,
	sample2,
	sample3,
	sample4,
	sample5,
	sample6
}