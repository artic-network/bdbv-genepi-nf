#!/bin/bash
# Print how many CPUs this process may actually use, capped at the number requested (first argument).
#
# Run inside the task, so it sees the real limits of a container: CPUs the container is pinned to
# (--cpuset-cpus, seen by nproc), a CPU quota (--cpus, in the cgroup), or the CPUs a Docker Desktop
# machine has been given. Falls back to the machine's CPUs outside a container.
requested="${1:-1}"

n=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)

quota=""
if [ -r /sys/fs/cgroup/cpu.max ]; then                      # cgroup v2
    read -r q p < /sys/fs/cgroup/cpu.max
    if [ "$q" != "max" ] && [ -n "$p" ] && [ "$p" -gt 0 ] 2>/dev/null; then
        quota=$(( q / p ))
    fi
elif [ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us ]; then       # cgroup v1
    q=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us)
    p=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)
    if [ "$q" -gt 0 ] 2>/dev/null && [ "$p" -gt 0 ] 2>/dev/null; then
        quota=$(( q / p ))
    fi
fi
if [ -n "$quota" ] && [ "$quota" -lt "$n" ]; then
    n=$quota
fi

[ "$n" -gt "$requested" ] && n=$requested
[ "$n" -lt 1 ] && n=1
echo "$n"
