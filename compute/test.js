let orbit = require('../index.js');
let mapOrbit = orbit.mapOrbit;

let testOrbit = require('./test_orbit.js');
let esc_sample = require('./test_escobal.js');

let gauss = require('./gauss.js');

let mu =3.9725202731797994e+14; //как у Эскобала


function different(res, sam){
	//console.log(sam);
	//console.log(res);
	let keys = ["Om", "om", "i", "p", "eps"];
	
	return keys.reduce((akk, key)=>(akk[key]=(res[key]-sam[key])/sam[key],akk), {});
}

function control(dif, top){
	for(let key of dif){
		if(dif[key]>top) return false;
	}
	return true;
}

function healtEsc(sam){
	sam.p = sam.a*(1-sam.eps**2);
	return sam;
}

function test(sam){
	let v1 = gauss(mu, sam.r1, sam.r2, sam.tau, 1, 0.0001)
	console.log(sam.r1, v1);
	let res = mapOrbit(mu, sam.r1, v1);

	let dif = different(res, sam);
	console.log(dif);
}


//test(esc_sample.sample1);
//[...testOrbit.tasks(mu, [esc_sample.sample1], Math.PI/8)];
let df = Math.PI/8, f1 = -Math.PI;
let orb = healtEsc(esc_sample.sample1);
let task = {...orb, ...testOrbit.twoR(mu, healtEsc(esc_sample.sample1), f1, f1+df)};


//test(healtEsc(esc_sample.sample1));
test(task);
