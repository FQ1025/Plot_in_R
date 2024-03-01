#!/bin/awk -f
BEGIN{
OFS="\t"
}
{
print $1,$2,$3,$4":"$5"-"$6","$7
}
