
dos2unix ./build_linux_distrib.sh

bash -i "build_linux_distrib.sh"

pause

rem rpmbuild -ba dsr-linux-1.2.11.spec