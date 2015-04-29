mata:
mata set matastrict on
	void function fe_parse(string scalar absvar) {
	
		

	
	// Preparse absvars because tokenget is stupid and splits foreign#c.(weight length)
	// Another tokenget bug is when dealing with pchars; it strips the spaces, so we need another workaround


	t = tokeninit("", "", ("()"))
	tokenset(t, absvar_names)
	absvars = tokengetall(t)
	for (g=1; g<=cols(absvars); g++) {
		absvars[g] = invtokens(tokens(absvars[g],"@"), ",")
	}
	absvars' 
	strlen(absvars')
	absvar_names
	invtokens(absvars)
	invtokens(absvars, "")
	invtokens(absvars, "&")
	123456
	


	t = tokeninit(" ", ("=", "#", "##", "i.", "c."), ("()") ) 
	t = tokeninit(" ", "", ("()"))
	// whitespace chars, parsing chars (whitespace that are also content), quote chars
	G = cols(absvars)
	//absvar_names
	//absvars
	G
	for (g=1; g<=G; g++) {
		fe_parse(S.fixed_effects[g], absvars[g])
	}

	S.foobar = 2
	S.spam = 7
	}
end
