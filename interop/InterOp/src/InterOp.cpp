#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <Rcpp.h>
using namespace Rcpp;

/***************************************
 *
 * read extraction metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::DataFrame readExtractionMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number (3)
	BYTE length;	// byte 1: length of each record
	// bytes (N * 38 + 2) - (N * 38 + 39): record: (N is the record index)
	#pragma pack(push, 1)
	struct ExtractionMetrics {
		uint16_t lane;	// 2 bytes: lane number
		uint16_t tile;	// 2 bytes: tile number
		uint16_t cycle;	// 2 bytes: cycle number
		//unsigned char padding[2];		// here the struct should be padded, but the file is not
		float    fwhmA;	// 4 bytes: fwhm scores for channel A respectively
		float    fwhmC;	// 4 bytes: fwhm scores for channel C respectively
		float    fwhmG;	// 4 bytes: fwhm scores for channel G respectively
		float    fwhmT;	// 4 bytes: fwhm scores for channel T respectively
		uint16_t intA;	// 2 bytes: intensities for channel A respectively
		uint16_t intC;	// 2 bytes: intensities for channel C respectively
		uint16_t intG;	// 2 bytes: intensities for channel G respectively
		uint16_t intT;	// 2 bytes: intensities for channel T respectively
		uint64_t datetime;	// 8 bytes: date/time of CIF creation (.Net timestamp (100 nanosec tics from 01-01-0001))
		};
	#pragma pack(pop)
	ExtractionMetrics reg;
	
	/*
	 * output data structures: vectors that will be put together into a df
	 */
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 2) / sizeof(reg) + 1;

	Rcpp::IntegerVector lane(l);
	Rcpp::IntegerVector tile(l);
	Rcpp::IntegerVector cycle(l);
	Rcpp::NumericVector fwhmA(l);
	Rcpp::NumericVector fwhmC(l);
	Rcpp::NumericVector fwhmG(l);
	Rcpp::NumericVector fwhmT(l);
	Rcpp::IntegerVector intA(l);
	Rcpp::IntegerVector intC(l);
	Rcpp::IntegerVector intG(l);
	Rcpp::IntegerVector intT(l);
	Rcpp::NumericVector datetime(l);	// 64 bits number. Unset the 2 most significant bits to get the 100 nanosec ticks since 01-01-0001

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);
	r = fread(&length, 1, 1, file);

	// loop over the registers
	int i = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		lane[i]      = reg.lane;
		tile[i]      = reg.tile;
		cycle[i]     = reg.cycle;
		fwhmA[i]     = reg.fwhmA;
		fwhmC[i]     = reg.fwhmC;
		fwhmG[i]     = reg.fwhmG;
		fwhmT[i]     = reg.fwhmT;
		intA[i]      = reg.intA;
		intC[i]      = reg.intC;
		intG[i]      = reg.intG;
		intT[i]      = reg.intT;
		datetime[i]  = reg.datetime & 0x3FFFFFFFFFFFFFFF; // remove the first 2 bits of this 64bit number (useless flags)
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane")     = lane,
		Rcpp::Named("tile")     = tile,
		Rcpp::Named("cycle")    = cycle,
		Rcpp::Named("fwhmA")    = fwhmA,
		Rcpp::Named("fwhmC")    = fwhmC,
		Rcpp::Named("fwhmG")    = fwhmG,
		Rcpp::Named("fwhmT")    = fwhmT,
		Rcpp::Named("intA")     = intA,
		Rcpp::Named("intC")     = intC,
		Rcpp::Named("intG")     = intG,
		Rcpp::Named("intT")     = intT,
		Rcpp::Named("datetime") = datetime);

	return df;
}

/***************************************
 *
 * read quality metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::List readQualityMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number (3)
	BYTE length;	// byte 1: length of each record
	// bytes (N * 206 + 2) - (N * 206 + 207): record: (N is the record index)
	#pragma pack(push, 1)
	struct QualityMetrics {
		uint16_t lane;	// 2 bytes: lane number
		uint16_t tile;	// 2 bytes: tile number
		uint16_t cycle;	// 2 bytes: metric cycle
		uint32_t nclust[50];	// number of clusters assigned score Q1 through Q50
		};
 	#pragma pack(pop)
 	QualityMetrics reg;
	
	/*
	 * output data structures: vectors that will be put together into a df
	 */
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 2) / sizeof(reg) + 1;

	Rcpp::IntegerVector lane(l);
	Rcpp::IntegerVector tile(l);
	Rcpp::IntegerVector cycle(l);
	Rcpp::IntegerMatrix nclust(l,50);

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);
	r = fread(&length, 1, 1, file);

	// loop over the registers
	int i,j = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		lane[i]  = reg.lane;
		tile[i]  = reg.tile;
		cycle[i] = reg.cycle;
		for(j=0; j<50; j++) {
			nclust(i,j) = reg.nclust[j];
		}
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane")  = lane,
		Rcpp::Named("tile")  = tile,
		Rcpp::Named("cycle") = cycle);

	Rcpp::List li = Rcpp::List::create(
		Rcpp::Named("key")   = df,
		Rcpp::Named("nclust")= nclust);

	return li;
}

/***************************************
 *
 * read error metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::DataFrame readErrorMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number (3)
	BYTE length;	// byte 1: length of each record
	// bytes (N * 30 + 2) - (N * 30 + 31): record: (N is the record index)
	#pragma pack(push, 1)
	struct ErrorMetrics {
		uint16_t lane;	// 2 bytes: lane number
		uint16_t tile;	// 2 bytes: tile number
		uint16_t cycle;	// 2 bytes: cycle number
		//unsigned char padding[2];		// here the struct should be padded, but the file is not
		float    erate;	// 4 bytes: error rate
		uint32_t n;		// 4 bytes: number of perfect reads
		uint32_t n1e;	// 4 bytes: number of reads with 1 error
		uint32_t n2e;	// 4 bytes: number of reads with 2 errors
		uint32_t n3e;	// 4 bytes: number of reads with 3 errors
		uint32_t n4e;	// 4 bytes: number of reads with 4 errors
		};
	#pragma pack(pop)
	ErrorMetrics reg;
	
	/*
	 * output data structures: vectors that will be put together into a df
	*/
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 2) / sizeof(reg) + 1;

	Rcpp::IntegerVector lane(l);
	Rcpp::IntegerVector tile(l);
	Rcpp::IntegerVector cycle(l);
	Rcpp::NumericVector erate(l);
	Rcpp::IntegerVector n(l);
	Rcpp::IntegerVector n1e(l);
	Rcpp::IntegerVector n2e(l);
	Rcpp::IntegerVector n3e(l);
	Rcpp::IntegerVector n4e(l);

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);
	r = fread(&length, 1, 1, file);

	// loop over the registers
	int i = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		lane[i]  = reg.lane;
		tile[i]  = reg.tile;
		cycle[i] = reg.cycle;
		erate[i] = reg.erate;
		n[i]     = reg.n;
		n1e[i]   = reg.n1e;
		n2e[i]   = reg.n2e;
		n3e[i]   = reg.n3e;
		n4e[i]   = reg.n4e;
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane") = lane,
	    Rcpp::Named("tile") = tile,
	    Rcpp::Named("cycle") = cycle,
	    Rcpp::Named("erate") = erate,
	    Rcpp::Named("n") = n,
	    Rcpp::Named("n1e") = n1e,
	    Rcpp::Named("n2e") = n2e,
	    Rcpp::Named("n3e") = n3e,
	    Rcpp::Named("n4e") = n4e);

	return df;
}

/***************************************
 *
 * read tile metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::DataFrame readTileMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number (3)
	BYTE length;	// byte 1: length of each record
	// bytes (N * 10 + 2) - (N * 10 + 11): record: (N is the record index)
	#pragma pack(push, 1)
	struct TileMetrics {
		uint16_t lane;	// 2 bytes: lane number
		uint16_t tile;	// 2 bytes: tile number
		uint16_t code;	// 2 bytes: metric code
		float    value;	// 4 bytes: metric value
		};
	/* possible metric codes are:
	 * code 100: cluster density (k/mm2)
	 * code 101: cluster density passing filters (k/mm2)
	 * code 102: number of clusters
	 * code 103: number of clusters passing filters
	 * code (200 + (N – 1) * 2): phasing for read N
	 * code (201 + (N – 1) * 2): prephasing for read N
	 * code (300 + N – 1): percent aligned for read N
	 * code 400: control lane */
 	#pragma pack(pop)
 	TileMetrics reg;
	
	/*
	 * output data structures: vectors that will be put together into a df
	 */
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 2) / sizeof(reg) + 1;

	Rcpp::IntegerVector lane(l);
	Rcpp::IntegerVector tile(l);
	Rcpp::IntegerVector code(l);
	Rcpp::NumericVector value(l);

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);
	r = fread(&length, 1, 1, file);

	// loop over the registers
	int i = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		lane[i]  = reg.lane;
		tile[i]  = reg.tile;
		code[i]  = reg.code;
		value[i] = reg.value;
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane") = lane,
		Rcpp::Named("tile") = tile,
		Rcpp::Named("code") = code,
		Rcpp::Named("value")= value);

	return df;
}

/***************************************
 *
 * read corrected intensity metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::DataFrame readCorrectedIntMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number (3)
	BYTE length;	// byte 1: length of each record
	// bytes (N * 48 + 2) - (N * 48 + 49): record: (N is the record index)
	#pragma pack(push, 1)
	struct CorrectedIntMetrics {
		uint16_t lane;			// 2 bytes: lane number
		uint16_t tile;      	// 2 bytes: tile number
		uint16_t cycle;     	// 2 bytes: cycle number
		uint16_t avgint;    	// 2 bytes: average intensity
		uint16_t avgintA;   	// 2 bytes: average corrected int for channel A
		uint16_t avgintC;   	// 2 bytes: average corrected int for channel C
		uint16_t avgintG;   	// 2 bytes: average corrected int for channel G
		uint16_t avgintT;   	// 2 bytes: average corrected int for channel T
		uint16_t avgintclA; 	// 2 bytes: average corrected int for called clusters for base A
		uint16_t avgintclC; 	// 2 bytes: average corrected int for called clusters for base C
		uint16_t avgintclG; 	// 2 bytes: average corrected int for called clusters for base G
		uint16_t avgintclT; 	// 2 bytes: average corrected int for called clusters for base T
		float    bcNC;      	// 4 bytes: number of base calls for No Call
		float    bcA;       	// 4 bytes: number of base calls for channel A
		float    bcC;       	// 4 bytes: number of base calls for channel C
		float    bcG;       	// 4 bytes: number of base calls for channel G
		float    bcT;       	// 4 bytes: number of base calls for channel T
		float    srratio;   	// 4 bytes: signal to noise ratio
		};
	#pragma pack(pop)
	CorrectedIntMetrics reg;
	
	/*
	 * output data structures: vectors that will be put together into a df
	 */
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 2) / sizeof(reg) + 1;

	Rcpp::IntegerVector lane(l);	   
	Rcpp::IntegerVector tile(l);     
	Rcpp::IntegerVector cycle(l);    
	Rcpp::IntegerVector avgint(l);   
	Rcpp::IntegerVector avgintA(l);  
	Rcpp::IntegerVector avgintC(l);  
	Rcpp::IntegerVector avgintG(l);  
	Rcpp::IntegerVector avgintT(l);  
	Rcpp::IntegerVector avgintclA(l);
	Rcpp::IntegerVector avgintclC(l);
	Rcpp::IntegerVector avgintclG(l);
	Rcpp::IntegerVector avgintclT(l);
	Rcpp::NumericVector bcNC(l);     
	Rcpp::NumericVector bcA(l);      
	Rcpp::NumericVector bcC(l);      
	Rcpp::NumericVector bcG(l);      
	Rcpp::NumericVector bcT(l);      
	Rcpp::NumericVector srratio(l);  

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);
	r = fread(&length, 1, 1, file);

	// loop over the registers
	int i = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		lane[i]     = reg.lane;
		tile[i]     = reg.tile;
		cycle[i]    = reg.cycle;
		avgint[i]   = reg.avgint;
		avgintA[i]  = reg.avgintA;
		avgintC[i]  = reg.avgintC;
		avgintG[i]  = reg.avgintG;
		avgintT[i]  = reg.avgintT;
		avgintclA[i]= reg.avgintclA;
		avgintclC[i]= reg.avgintclC;
		avgintclG[i]= reg.avgintclG;
		avgintclT[i]= reg.avgintclT;
		bcNC[i]     = reg.bcNC;
		bcA[i]      = reg.bcA;
		bcC[i]      = reg.bcC;
		bcG[i]      = reg.bcG;
		bcT[i]      = reg.bcT;
		srratio[i]  = reg.srratio;
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane")      = lane,
		Rcpp::Named("tile")      = tile,
		Rcpp::Named("cycle")     = cycle,
		Rcpp::Named("avgint")    = avgint,
		Rcpp::Named("avgintA")   = avgintA,
		Rcpp::Named("avgintC")   = avgintC,
		Rcpp::Named("avgintG")   = avgintG,
		Rcpp::Named("avgintT")   = avgintT,
		Rcpp::Named("avgintclA") = avgintclA,
		Rcpp::Named("avgintclC") = avgintclC,
		Rcpp::Named("avgintclG") = avgintclG,
		Rcpp::Named("avgintclT") = avgintclT,
		Rcpp::Named("bcNC")      = bcNC,
		Rcpp::Named("bcA")       = bcA,
		Rcpp::Named("bcC")       = bcC,
		Rcpp::Named("bcG")       = bcG,
		Rcpp::Named("bcT")       = bcT,
		Rcpp::Named("srratio")   = srratio);

	return df;
}

/***************************************
 *
 * read control metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::DataFrame readControlMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number
	// the fixed length part of the register
	#pragma pack(push, 1)
	struct ControlMetrics {
		uint16_t lane;	// 2 bytes: lane number
		uint16_t tile;	// 2 bytes: tile number
		uint16_t read;	// 2 bytes: read number
	};
	#pragma pack(pop)
	ControlMetrics reg;
	// the variable part of the register
	uint16_t reg_control_l;	// 2 bytes: number bytes X for control name
	char     reg_control[256];	//X bytes: control name string (string in UTF8Encoding)
	uint16_t reg_index_l;	// 2 bytes: number bytes Y for index name
	char     reg_index[256];// Y bytes: index name string (string in UTF8Encoding)
	uint32_t reg_nclust;	// 4 bytes: # clusters identified as control
	
	/*
	 * output data structures: vectors that will be put together into a df
	 */
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 1) / (sizeof(reg) + 10) + 1;	// 10 is the minimum size of the variable part

	Rcpp::IntegerVector   lane(l);
	Rcpp::IntegerVector   tile(l);
	Rcpp::IntegerVector   read(l);
	Rcpp::CharacterVector control(l);
	Rcpp::CharacterVector index(l);
	Rcpp::IntegerVector   nclust(l);

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);

	// loop over the registers
	std::string s;	// for char* to string conversion
	int i = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		r = fread(&reg_control_l, sizeof(reg_control_l), 1, file);
		r = fread(&reg_control, reg_control_l, 1, file);
		r = fread(&reg_index_l, sizeof(reg_index_l), 1, file);
		r = fread(&reg_index, reg_index_l, 1, file);
		r = fread(&reg_nclust, sizeof(reg_nclust), 1, file);
		lane[i]   = reg.lane;
		tile[i]   = reg.tile;
		read[i]   = reg.read;
		reg_control[reg_control_l] = '\0';	// put the end of string mark
		s         = reg_control;	// implicit char* to string conversion
		control[i]= s;
		reg_index[reg_index_l] = '\0';	// put the end of string mark
		s         = reg_index;		// implicit char* to string conversion
		index[i]  = s;
		nclust[i] = reg_nclust;
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	// subset only the 'i' read registers. Special subseting for CharacterVector
	Rcpp::IntegerVector   lane_i    = lane[seq(0,i-1)]; 
	Rcpp::IntegerVector   tile_i    = tile[seq(0,i-1)]; 
	Rcpp::IntegerVector   read_i    = read[seq(0,i-1)]; 
	Rcpp::CharacterVector control_i(i);
	std::copy(control.begin(), control.begin() + i, control_i.begin()); 
	Rcpp::CharacterVector index_i(i);
	std::copy(index.begin()  , index.begin() + i  , index_i.begin()); 
	Rcpp::IntegerVector   nclust_i  = nclust[seq(0,i-1)]; 

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane")   = lane_i,
		Rcpp::Named("tile")   = tile_i,
		Rcpp::Named("read")   = read_i,
		Rcpp::Named("control")= control_i,
		Rcpp::Named("index")  = index_i,
		Rcpp::Named("nclust") = nclust_i);

	return df;
}

/***************************************
 *
 * read image metrics
 *
 ***************************************/
// [[Rcpp::export]]
Rcpp::DataFrame readImageMetrics(CharacterVector f) {
	
	std::string fx = as<std::string>(f[0]);
	
	/*
	 * register definition
	 */
	typedef unsigned char BYTE;
	BYTE version;	// byte 0: file version number (1)
	BYTE length;	// byte 1: length of each record
	// bytes (N * 12 + 2) - (N * 12 + 13): record: (N is the record index)
	#pragma pack(push, 1)
	struct ImageMetrics {
		uint16_t lane;	// 2 bytes: lane number
		uint16_t tile;	// 2 bytes: tile number
		uint16_t cycle;	// 2 bytes: cycle number
		uint16_t channelid;	// 2 bytes: channel id where 0=A, 1=C, 2=G, 3=T
		uint16_t mincont;	// 2 bytes: min contrast value for image
		uint16_t maxcont;	// 2 bytes: max contrast value for image
		};
	#pragma pack(pop)
	ImageMetrics reg;
	
	/*
	 * output data structures: vectors that will be put together into a df
	*/
	struct stat st;
	if(stat(fx.c_str(), &st) != 0) {
		stop("Could not open specified file");
	}
	int l = (st.st_size - 2) / sizeof(reg) + 1;

	Rcpp::IntegerVector lane(l);
	Rcpp::IntegerVector tile(l);
	Rcpp::IntegerVector cycle(l);
	Rcpp::IntegerVector channelid(l);
	Rcpp::IntegerVector mincont(l);
	Rcpp::IntegerVector maxcont(l);

	/*
	 * read input and feed the output structure
	 */
	// open input file handle
	FILE *file = NULL;
	if ((file = fopen(fx.c_str(), "rb")) == NULL) {
		stop("Could not open specified file");
	}

	// read version and register length
	int r;
	r = fread(&version, 1, 1, file);
	r = fread(&length, 1, 1, file);

	// loop over the registers
	int i = 0;
	while(!feof(file)) {
		r = fread(&reg, sizeof(reg), 1, file);
		lane[i]     = reg.lane;
		tile[i]     = reg.tile;
		cycle[i]    = reg.cycle;
		channelid[i]= reg.channelid;
		mincont[i]  = reg.mincont;
		maxcont[i]  = reg.maxcont;
		i++;
	}
	
	// close the file and return the data frame
	fclose(file);

	Rcpp::DataFrame df = Rcpp::DataFrame::create(
		Rcpp::Named("lane")     = lane,
		Rcpp::Named("tile")     = tile,
		Rcpp::Named("cycle")    = cycle,
		Rcpp::Named("channelid")= channelid,
		Rcpp::Named("mincont")  = mincont,
		Rcpp::Named("maxcont")  = maxcont);

	return df;
}
