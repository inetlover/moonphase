<?php

require('Classes/Moonphase.php');

$moonphase = new \TTREE\moonphase();

// phaseHunting() Example
print "<pre>";
print "Example: phaseHunting()\n";
print "-----------------------\n";
$moonphase->calculatePhaseHunting();
print "\n\n";


// findPhaseList() Example
print "Example: findPhaseList()\n";
print "------------------------\n";
$startDate = new DateTime("2012-01-01 00:00:00 CEST");
$endDate = new DateTime("2012-12-01 00:00:00 CEST");
$moonphase->calculatePhaseList($startDate, $endDate);
print "\n\n";


// calculatePhase() Example
$date = new DateTime('2012-01-01 00:00:00 CEST');

print "Example: calculatePhase()\n";
print "-------------------------\n";
$moonphase->calculatePhaseByDate($date);
print "</pre>\n\n";


?>
