let {sign, sqrt} = Math;
function common(r1, r2){
	let U1=r1.ort(),
		U2=r2.ort(),
		dfi_cos=U1.smul(U2),
		W = r1.vmul(r2).ort(), //*s
		dfi_sin = sign(W.z)*sqrt(1-dfi_cos**2);
	return {
		U1, U2, dfi_cos, dfi_sin, W
	}
}

module.exports = common;