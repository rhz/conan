/*
 * Copyright (C) 2003  Columbia Genome Center
 * All Rights Reserved.
 *
 * parseargs -- parses arguments from command line
 *
 * $Id: parseargs.h,v 1.1 2005/03/22 17:41:33 kw2110 Exp $
 */
#define	ARGBEGIN \
	for (argc--, argv++; \
	    *argv != NULL && argv[0][0] == '-' && argv[0][1] != '\0'; \
	    argc--, argv++) \
	{ \
		int _argc; \
		char *_tmp, *_arg = argv[0] + 1; \
		if (_arg[0] == '-' && _arg[1] == '\0') { \
			argc--; \
			argv++; \
			break; \
		} \
		_tmp = 0; \
		while ((_argc = *_arg++) != '\0') \
			switch (_argc)
#define	ARGC()	*argv
#define	ARGF()	(_tmp = _arg, _arg = "", \
		    (_tmp[0] != '\0') ? \
			_tmp : \
			(argv[1] != NULL) ? (argc--, *++argv) : "")
#define	ARGEND	}
