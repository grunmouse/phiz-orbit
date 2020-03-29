module.exports = function newton_raffson(F, dF, x, eps){
	let del;
	do{
		del = F(x)/dF(x);
		x = x - del;
	}while(Math.abs(del)>eps);
	return x;
}