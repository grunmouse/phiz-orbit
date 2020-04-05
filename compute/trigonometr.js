
function sin2a(sina, cosa){
	return 2*sina*cosa;
}

function cos2a(sina, cosa){
	if(isNaN(sina)){
		return 2*cosa**2-1;
	}
	else if(isNaN(cosa)){
		return 1-2*sina**2
	}
}

function cosa2(cosa){
	return sqrt((1+cosa)/2)
}


module.exports = {
	sin2a,
	cos2a,
	cosa2
}