<?php

require('moonphase.inc.php');

$moonphase = new \TTREE\moonphase();

// phasehunt() Example
print "<pre>";
print "Example: phasehunt()\n";
$moonphase->do_phasehunt();
print "\n\n";


// phaselist() Example
print "Example: phaselist()\n";
$start = strtotime("2012-01-01 00:00:00 CEST");
$stop = strtotime("2012-12-31 00:00:00 CEST");
$moonphase->do_phaselist($start, $stop);
print "\n\n";


// phase() Example
$date = "2012-01-01";
$time = "00:00:00";
$tzone = "CEST";
print "Example: phase() ($date $time $tzone)\n";
$moonphase->do_phase($date, $time, $tzone);
print "</pre>\n\n";


?>
