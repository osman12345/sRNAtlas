#!/bin/bash
# Launch script for sRNA-seq Analysis WebTool

# Set environment
export PYTHONPATH="${PYTHONPATH}:$(dirname "$0")"

# Check if streamlit is installed
if ! command -v streamlit &> /dev/null; then
    echo "Streamlit not found. Installing..."
    pip install streamlit streamlit-option-menu
fi

# Default port
PORT=${1:-8501}

echo "=============================================="
echo "  sRNA-seq Analysis WebTool"
echo "=============================================="
echo ""
echo "Starting server on port $PORT..."
echo "Open http://localhost:$PORT in your browser"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Run streamlit
streamlit run app/main.py --server.port $PORT --server.headless true
