#!/bin/bash

# Script to update snapshots for all nf-test test files

echo "Updating snapshots for all tests..."
echo

test_files=("default" "test_donor_match" "test_genetic" "test_hashing")

for test_file in "${test_files[@]}"; do
    if [ "$test_file" != "default" ]; then
        test_profile="$test_file"
    else
        test_profile="test"
    fi

    command="nf-test test tests/${test_file}.nf.test --profile ${test_profile},docker --update-snapshot"

    echo "Updating snapshot for: $test_file"
    echo "Running: ${command}"

    eval "$command"

    echo "✓ Snapshot updated for: $test_file"
    echo
done

echo "All snapshots have been updated!"
