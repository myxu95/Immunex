#!/bin/bash
# Monitor batch allostery analysis progress

LOG_FILE="/home/xumy/work/development/AfterMD/batch_allostery_8core.log"
PROCESS_ID=1869918

echo "========================================================================"
echo "Batch Allostery Analysis - 8 Core Parallel Processing Monitor"
echo "========================================================================"
echo ""

# Check if process is running
if ps -p $PROCESS_ID > /dev/null 2>&1; then
    echo "✓ Process is RUNNING (PID: $PROCESS_ID)"

    # Get process runtime
    RUNTIME=$(ps -p $PROCESS_ID -o etime= | xargs)
    echo "✓ Runtime: $RUNTIME"

    # Get CPU usage
    CPU=$(ps -p $PROCESS_ID -o %cpu= | xargs)
    echo "✓ CPU usage: ${CPU}%"

    # Get memory usage
    MEM=$(ps -p $PROCESS_ID -o %mem= | xargs)
    echo "✓ Memory usage: ${MEM}%"
else
    echo "✗ Process is NOT running (may have completed or failed)"
fi

echo ""
echo "------------------------------------------------------------------------"
echo "Progress Information:"
echo "------------------------------------------------------------------------"

# Extract progress from log
PROGRESS_LINE=$(grep -a "Progress:" "$LOG_FILE" | tail -1)
if [ -n "$PROGRESS_LINE" ]; then
    echo "$PROGRESS_LINE"
else
    echo "No progress information available yet..."
fi

echo ""
echo "------------------------------------------------------------------------"
echo "Recently Completed Tasks (last 5):"
echo "------------------------------------------------------------------------"
grep -a "Completed" "$LOG_FILE" | tail -5

echo ""
echo "------------------------------------------------------------------------"
echo "Currently Processing (last 10 frame updates):"
echo "------------------------------------------------------------------------"
grep -a "Frame.*/" "$LOG_FILE" | tail -10

echo ""
echo "------------------------------------------------------------------------"
echo "Error Check:"
echo "------------------------------------------------------------------------"
ERROR_COUNT=$(grep -ac "ERROR\|Error\|Failed" "$LOG_FILE" 2>/dev/null || echo "0")
echo "Total errors found: $ERROR_COUNT"

if [ "$ERROR_COUNT" -gt 0 ]; then
    echo ""
    echo "Recent errors:"
    grep -a "ERROR\|Error\|Failed" "$LOG_FILE" | tail -5
fi

echo ""
echo "========================================================================"
echo "Commands:"
echo "  Watch live: tail -f $LOG_FILE"
echo "  Kill process: kill $PROCESS_ID"
echo "  Check outputs: ls -lh output/allostery_analysis/"
echo "========================================================================"
