/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.limiters;

import std.math : abs, fmax, fmin, sgn, sqrt;
import std.meta : aliasSeqOf;

alias limiterList = aliasSeqOf!(["firstOrderUpwindS", "laxWendroffS", "minmodS", "harmonicS", "geometricS", "superbeeS"]);

double firstOrderUpwindS(double Sp, double Sm)
{
	return 0.0;
}

double laxWendroffS(double Sp, double Sm)
{
	return Sp;
}

double minmodS(double Sp, double Sm)
{
	return fmax(0, fmin(abs(Sm), abs(Sp)))*sgn(Sm + Sp);
}

double harmonicS(double Sp, double Sm)
{
	return (Sm*abs(Sp) + abs(Sm)*Sp)/(abs(Sm) + abs(Sp) + 10e-12);
}

double geometricS(double Sp, double Sm)
{
	return sqrt( (Sm*Sp + abs(Sm*Sp))/2)*sgn(Sm + Sp);
}

double superbeeS(double Sp, double Sm)
{
	return fmax(0, fmax(fmin(2*abs(Sm), abs(Sp)), fmin(abs(Sm), 2*abs(Sp))))*sgn(Sm + Sp);
}