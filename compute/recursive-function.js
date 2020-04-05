function recursiveFunction(store, fun){
	return function me(n){
		if(n<0){
			throw new Error('Invalid argument '+ n);
		}
		if(n<store.length){
			return store[n];
		}
		else{
			return store[n] = fun(me, n);
		}
	}
}

module.exports = recursiveFunction;