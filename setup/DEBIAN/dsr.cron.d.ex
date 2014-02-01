#
# Regular cron jobs for the dsr package
#
0 4	* * *	root	[ -x /usr/bin/dsr_maintenance ] && /usr/bin/dsr_maintenance
