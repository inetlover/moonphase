<?php
namespace TTREE;

class moonphase {

	// Error definition
	const ERR_UNDEF = -1;

	// Astronomical constants
	const EPOCH = 2444238.5;

	// Constants defining the Sun's apparent orbit.
	const ELONGE = 278.833540; // ecliptic longitude of the Sun at epoch 1980.0
	const ELONGP = 282.596403; // ecliptic longitude of the Sun at perigee
	const ECCENT = 0.016718; // eccentricity of Earth's orbit
	const SUNSMAX = 1.495985e8; // semi-major axis of Earth's orbit, km
	const SUNANGSIZ = 0.533128; // sun's angular size, degrees, at semi-major axis distance

	// Elements of the Moon's orbit, epoch 1980.0.
	const MMLONG = 64.975464; // moon's mean longitude at the epoch
	const MMLONGP = 64.975464; // ecliptic longitude of the Sun at perigee
	const MLNODE = 151.950429; // mean longitude of the node at the epoch
	const MINC = 151.950429; // inclination of the Moon's orbit
	const MECC = 5.145396; // inclination of the Moon's orbit
	const MANGSIZ = 0.054900; // eccentricity of the Moon's orbit
	const MSMAX = 384401.0; // semi-major axis of Moon's orbit in km
	const MPARALLAX = 0.9507; // parallax at distance a from Earth
	const SYNMONTH = 29.53058868; // synodic month (new Moon to new Moon)


	/**
	 * Phase Hunt
	 */
	function do_phasehunt() {
		$phases = $this->phasehunt();

		print date("D M j G:i:s T Y", $phases[0]) . "\n";
		print date("D M j G:i:s T Y", $phases[1]) . "\n";
		print date("D M j G:i:s T Y", $phases[2]) . "\n";
		print date("D M j G:i:s T Y", $phases[3]) . "\n";
		print date("D M j G:i:s T Y", $phases[4]) . "\n";
	}


	/**
	 * Phase list
	 *
	 * @param $start
	 * @param $stop
	 */
	function do_phaselist($start, $stop) {
		$name = array("New Moon", "First quarter", "Full moon", "Last quarter");
		$times = $this->phaselist($start, $stop);

		foreach ($times as $time) {

			// First element is the starting phase (see $name).
			if ($time == $times[0]) {
				print $name["$times[0]"] . "\n";
			}
			else {
				print date("D M j G:i:s T Y", $time) . "\n";
			}
		}
	}

	/**
	 * Calculate phase
	 *
	 * @param int $date
	 * @param int $time
	 * @param string $tzone
	 */
	function do_phase($date, $time, $tzone) {
		$moondata = $this->phase(strtotime($date . ' ' . $time . ' ' . $tzone));

		$MoonIllum = $moondata[1];
		$MoonAge = $moondata[2];

		$phase = 'Waxing';
		if ($MoonAge > self::SYNMONTH / 2) {
			$phase = 'Waning';
		}

		// Convert $MoonIllum to percent and round to whole percent.
		$MoonIllum = round($MoonIllum, 2);
		$MoonIllum *= 100;
		if ($MoonIllum == 0) {
			$phase = "New Moon";
		}
		if ($MoonIllum == 100) {
			$phase = "Full Moon";
		}

		print "Moon Phase: $phase\n";
		print "Percent Illuminated: $MoonIllum%\n";
	}

	/**
	 * Extract sign
	 *
	 * @param float $arg
	 * @return int
	 */
	function sgn($arg) {
		return (($arg < 0) ? -1 : ($arg > 0 ? 1 : 0));
	}

	/**
	 * Fix angle
	 *
	 * @param float $arg
	 * @return float
	 */
	function fixangle($arg) {
		return ($arg - 360.0 * (floor($arg / 360.0)));
	}

	/**
	 * Degree to radiant conversion
	 *
	 * @param float $arg
	 * @return float mixed
	 */
	function torad($arg) {
		return ($arg * (pi() / 180.0));
	}

	/**
	 * Radiant to degree conversion
	 *
	 * @param float $arg
	 * @return float mixed
	 */
	function todeg($arg) {
		return ($arg * (180.0 / pi()));
	}

	/**
	 * Sin from radian conversion
	 *
	 * @param float $arg
	 * @return float
	 */
	function dsin($arg) {
		return (sin($this->torad($arg)));
	}

	/**
	 * Cos from degree
	 *
	 * @param float $arg
	 * @return float
	 */
	function dcos($arg) {
		return (cos($this->torad($arg)));
	}

	/**
	 * Convert internal date and time to astronomical Julian
	 * time (i.e. Julian date plus day fraction)
	 *
	 * @param int $timestamp
	 * @return float
	 */
	function jtime($timestamp) {
		$julian = ($timestamp / 86400) + 2440587.5; // (seconds / (seconds per day)) + julian date of epoch
		return $julian;
	}

	/**
	 * Convert Julian date to a UNIX epoch
	 *
	 * @param int $jday
	 * @return int
	 */
	function jdaytosecs($jday = 0) {
		$stamp = ($jday - 2440587.5) * 86400; // (juliandate - jdate of unix epoch) * (seconds per julian day)
		return $stamp;
	}

	/**
	 * Convert Julian date to year, month, day, which are
	 * returned via integer pointers to integers
	 *
	 * @param float $td
	 * @param int $yy
	 * @param int $mm
	 * @param int $dd
	 */
	function jyear($td, &$yy, &$mm, &$dd) {
		$td += 0.5; // astronomical to civil.
		$z = floor($td);
		$f = $td - $z;

		if ($z < 2299161.0) {
			$a = $z;
		}
		else {
			$alpha = floor(($z - 1867216.25) / 36524.25);
			$a = $z + 1 + $alpha - floor($alpha / 4);
		}

		$b = $a + 1524;
		$c = floor(($b - 122.1) / 365.25);
		$d = floor(365.25 * $c);
		$e = floor(($b - $d) / 30.6001);

		$dd = $b - $d - floor(30.6001 * $e) + $f;
		$mm = $e < 14 ? $e - 1 : $e - 13;
		$yy = $mm > 2 ? $c - 4716 : $c - 4715;
	}


	/**
	 * Calculates time of the mean new Moon for a given
	 * base date. This argument K to this function is the
	 * precomputed synodic month index, given by:
	 *
	 * K = (year - 1900) * 12.3685
	 *
	 * where year is expressed as a year and fractional year.
	 *
	 * @param float $sdate
	 * @param float $k
	 * @return float
	 */
	function meanphase($sdate, $k) {

		// Time in Julian centuries from 1900 January 0.5
		$t = ($sdate - 2415020.0) / 36525;
		$t2 = $t * $t; // Square for frequent use
		$t3 = $t2 * $t; // Cube for frequent use

		$nt1 = 2415020.75933 + self::SYNMONTH * $k
			+ 0.0001178 * $t2
			- 0.000000155 * $t3
			+ 0.00033 * $this->dsin(166.56 + 132.87 * $t - 0.009173 * $t2);

		return ($nt1);
	}

	/**
	 * Given a K value used to determine the mean phase of the
	 * new moon, and a phase selector (0.0, 0.25, 0.5, 0.75),
	 * obtain the true, corrected phase time.
	 *
	 * @param float $k
	 * @param float $phase
	 * @return float
	 */
	function truephase($k, $phase) {
		$apcor = 0;

		$k += $phase; // add phase to new moon time
		$t = $k / 1236.85; // time in Julian centuries from 1900 January 0.5
		$t2 = $t * $t; // square for frequent use
		$t3 = $t2 * $t; // cube for frequent use

		// mean time of phase
		$pt = 2415020.75933
			+ self::SYNMONTH * $k
			+ 0.0001178 * $t2
			- 0.000000155 * $t3
			+ 0.00033 * $this->dsin(166.56 + 132.87 * $t - 0.009173 * $t2);

		// Sun's mean anomaly
		$m = 359.2242
			+ 29.10535608 * $k
			- 0.0000333 * $t2
			- 0.00000347 * $t3;

		// Moon's mean anomaly
		$mprime = 306.0253
			+ 385.81691806 * $k
			+ 0.0107306 * $t2
			+ 0.00001236 * $t3;

		// Moon's argument of latitude
		$f = 21.2964
			+ 390.67050646 * $k
			- 0.0016528 * $t2
			- 0.00000239 * $t3;

		if (($phase < 0.01) || (abs($phase - 0.5) < 0.01)) {
			// Corrections for New and Full Moon.
			$pt += (0.1734 - 0.000393 * $t) * $this->dsin($m)
				+ 0.0021 * $this->dsin(2 * $m)
				- 0.4068 * $this->dsin($mprime)
				+ 0.0161 * $this->dsin(2 * $mprime)
				- 0.0004 * $this->dsin(3 * $mprime)
				+ 0.0104 * $this->dsin(2 * $f)
				- 0.0051 * $this->dsin($m + $mprime)
				- 0.0074 * $this->dsin($m - $mprime)
				+ 0.0004 * $this->dsin(2 * $f + $m)
				- 0.0004 * $this->dsin(2 * $f - $m)
				- 0.0006 * $this->dsin(2 * $f + $mprime)
				+ 0.0010 * $this->dsin(2 * $f - $mprime)
				+ 0.0005 * $this->dsin($m + 2 * $mprime);
			$apcor = 1;
		}
		elseif ((abs($phase - 0.25) < 0.01 || (abs($phase - 0.75) < 0.01))) {
			$pt += (0.1721 - 0.0004 * $t) * $this->dsin($m)
				+ 0.0021 * $this->dsin(2 * $m)
				- 0.6280 * $this->dsin($mprime)
				+ 0.0089 * $this->dsin(2 * $mprime)
				- 0.0004 * $this->dsin(3 * $mprime)
				+ 0.0079 * $this->dsin(2 * $f)
				- 0.0119 * $this->dsin($m + $mprime)
				- 0.0047 * $this->dsin($m - $mprime)
				+ 0.0003 * $this->dsin(2 * $f + $m)
				- 0.0004 * $this->dsin(2 * $f - $m)
				- 0.0006 * $this->dsin(2 * $f + $mprime)
				+ 0.0021 * $this->dsin(2 * $f - $mprime)
				+ 0.0003 * $this->dsin($m + 2 * $mprime)
				+ 0.0004 * $this->dsin($m - 2 * $mprime)
				- 0.0003 * $this->dsin(2 * $m + $mprime);
			if ($phase < 0.5) {
				// First quarter correction.
				$pt += 0.0028 - 0.0004 * $this->dcos($m) + 0.0003 * $this->dcos($mprime);
			}
			else {
				// Last quarter correction.
				$pt += -0.0028 + 0.0004 * $this->dcos($m) - 0.0003 * $this->dcos($mprime);
			}
			$apcor = 1;
		}
		if (!$apcor) {
			print "truephase() called with invalid phase selector ($phase).\n";
			exit(self::ERR_UNDEF);
		}
		return ($pt);
	}

	/**
	 * Find time of phases of the moon which surround the current
	 * date.  Five phases are found, starting and ending with the
	 * new moons which bound the current lunation
	 *
	 * @param int $time
	 * @return array
	 */
	function phasehunt($time = -1) {

		if (empty($time) || $time == -1) {
			$time = time();
		}
		$sdate = $this->jtime($time);
		$adate = $sdate - 45;
		$this->jyear($adate, $yy, $mm, $dd);
		$k1 = floor(($yy + (($mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685);
		$adate = $nt1 = $this->meanphase($adate, $k1);

		$k2 = $k1;

		while (1) {
			$adate += self::SYNMONTH;
			$k2 = $k1 + 1;
			$nt2 = $this->meanphase($adate, $k2);
			if (($nt1 <= $sdate) && ($nt2 > $sdate)) {
				break;
			}
			$nt1 = $nt2;
			$k1 = $k2;
		}

		return array($this->jdaytosecs($this->truephase($k1, 0.0)),
			$this->jdaytosecs($this->truephase($k1, 0.25)),
			$this->jdaytosecs($this->truephase($k1, 0.5)),
			$this->jdaytosecs($this->truephase($k1, 0.75)),
			$this->jdaytosecs($this->truephase($k2, 0.0))
		);
	}

	/**
	 * Find time of phases of the moon between two dates.
	 * Times (in & out) are seconds_since_1970
	 *
	 * @param int $sdate
	 * @param int $edate
	 * @return array
	 */
	function phaselist($sdate, $edate) {
		if (empty($sdate) || empty($edate)) {
			return array();
		}

		$sdate = $this->jtime($sdate);
		$edate = $this->jtime($edate);

		$phases = array();
		$d = $k = $yy = $mm = 0;

		$this->jyear($sdate, $yy, $mm, $d);
		$k = floor(($yy + (($mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685) - 2;

		while (1) {
			++$k;
			foreach (array(0.0, 0.25, 0.5, 0.75) as $phase) {
				$d = $this->truephase($k, $phase);
				if ($d >= $edate) {
					return $phases;
				}
				if ($d >= $sdate) {
					if (empty($phases)) {
						array_push($phases, floor(4 * $phase));
					}
					array_push($phases, $this->jdaytosecs($d));
				}
			}
		}

		return array();
	}

	/**
	 * Solve the equation of Kepler
	 *
	 * @param float $m
	 * @param float $ecc
	 * @return float
	 */
	function kepler($m, $ecc) {
		$EPSILON = 1e-6;
		$m = $this->torad($m);
		$e = $m;
		do {
			$delta = $e - $ecc * sin($e) - $m;
			$e -= $delta / (1 - $ecc * cos($e));
		} while (abs($delta) > $EPSILON);
		return ($e);
	}

	/**
	 * Calculate phase of moon as a fraction:
	 *
	 * The argument is the time for which the phase is requested,
	 * expressed as a Julian date and fraction.  Returns the terminator
	 * phase angle as a percentage of a full circle (i.e., 0 to 1),
	 * and stores into pointer arguments the illuminated fraction of
	 * the Moon's disc, the Moon's age in days and fraction, the
	 * distance of the Moon from the centre of the Earth, and the
	 * angular diameter subtended by the Moon as seen by an observer
	 * at the centre of the Earth.
	 *
	 * @param int $time
	 * @return array
	 */
	function phase($time = 0) {
		if (empty($time) || $time == 0) {
			$time = time();
		}
		$pdate = $this->jtime($time);

		//	my ($Day, $N, $M, $Ec, $Lambdasun, $ml, $MM, $MN, $Ev, $Ae, $A3, $MmP,
		//	   $mEc, $A4, $lP, $V, $lPP, $NP, $y, $x, $Lambdamoon, $BetaM,
		//	   $MoonAge, $MoonPhase,
		//	   $MoonDist, $MoonDFrac, $MoonAng, $MoonPar,
		//	   $F, $SunDist, $SunAng,
		//	   $mpfrac);

		// Calculation of the Sun's position.
		$Day = $pdate - self::EPOCH; // date within epoch
		$N = $this->fixangle((360 / 365.2422) * $Day); // mean anomaly of the Sun
		$M = $this->fixangle($N + self::ELONGE - self::ELONGP); // convert from perigee co-ordinates
		//   to epoch 1980.0
		$Ec = $this->kepler($M, self::ECCENT); // solve equation of Kepler
		$Ec = sqrt((1 + self::ECCENT) / (1 - self::ECCENT)) * tan($Ec / 2);
		$Ec = 2 * $this->todeg(atan($Ec)); // true anomaly
		$Lambdasun = $this->fixangle($Ec + self::ELONGP); // Sun's geocentric ecliptic longitude
		# Orbital distance factor.
		$F = ((1 + self::ECCENT * cos($this->torad($Ec))) / (1 - self::ECCENT * self::ECCENT));
		$SunDist = self::SUNSMAX / $F; // distance to Sun in km
		$SunAng = $F * self::SUNANGSIZ; // Sun's angular size in degrees


		// Calculation of the Moon's position.

		// Moon's mean longitude.
		$ml = $this->fixangle(13.1763966 * $Day + self::MMLONG);

		// Moon's mean anomaly.
		$MM = $this->fixangle($ml - 0.1114041 * $Day - self::MMLONGP);

		// Moon's ascending node mean longitude.
		$MN = $this->fixangle(self::MLNODE - 0.0529539 * $Day);

		// Evection.
		$Ev = 1.2739 * sin($this->torad(2 * ($ml - $Lambdasun) - $MM));

		// Annual equation.
		$Ae = 0.1858 * sin($this->torad($M));

		// Correction term.
		$A3 = 0.37 * sin($this->torad($M));

		// Corrected anomaly.
		$MmP = $MM + $Ev - $Ae - $A3;

		// Correction for the equation of the centre.
		$mEc = 6.2886 * sin($this->torad($MmP));

		// Another correction term.
		$A4 = 0.214 * sin($this->torad(2 * $MmP));

		// Corrected longitude.
		$lP = $ml + $Ev + $mEc - $Ae + $A4;

		// Variation.
		$V = 0.6583 * sin($this->torad(2 * ($lP - $Lambdasun)));

		// True longitude.
		$lPP = $lP + $V;

		// Corrected longitude of the node.
		$NP = $MN - 0.16 * sin($this->torad($M));

		// Y inclination coordinate.
		$y = sin($this->torad($lPP - $NP)) * cos($this->torad(self::MINC));

		// X inclination coordinate.
		$x = cos($this->torad($lPP - $NP));

		// Ecliptic longitude.
		$Lambdamoon = $this->todeg(atan2($y, $x));
		$Lambdamoon += $NP;

		// Ecliptic latitude.
		$BetaM = $this->todeg(asin(sin($this->torad($lPP - $NP)) * sin($this->torad(self::MINC))));


		// Calculation of the phase of the Moon.

		// Age of the Moon in degrees.
		$MoonAge = $lPP - $Lambdasun;

		// Phase of the Moon.
		$MoonPhase = (1 - cos($this->torad($MoonAge))) / 2;

		// Calculate distance of moon from the centre of the Earth.
		$MoonDist = (self::MSMAX * (1 - self::MECC * self::MECC)) / (1 + self::MECC * cos($this->torad($MmP + $mEc)));

		// Calculate Moon's angular diameter.
		$MoonDFrac = $MoonDist / self::MSMAX;
		$MoonAng = self::MANGSIZ / $MoonDFrac;

		// Calculate Moon's parallax.
		$MoonPar = self::MPARALLAX / $MoonDFrac;

		$pphase = $MoonPhase;
		$mage = self::SYNMONTH * ($this->fixangle($MoonAge) / 360.0);
		$dist = $MoonDist;
		$angdia = $MoonAng;
		$sudist = $SunDist;
		$suangdia = $SunAng;
		$mpfrac = $this->fixangle($MoonAge) / 360.0;

		return array($mpfrac, $pphase, $mage, $dist, $angdia, $sudist, $suangdia);
	}
}

?>
