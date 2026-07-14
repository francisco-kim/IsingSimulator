// Canvas + animation-loop glue for the Ising simulator.
// Loaded both via JSHost.ImportAsync (for [JSImport] fast-path rendering) and
// via IJSRuntime module import (for the requestAnimationFrame loop); ES modules
// are singletons per URL, so both share this state.

const canvases = new Map();

export function initCanvas(id, width, height) {
  const canvas = document.getElementById(id);
  if (!canvas) {
    return;
  }
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d', { alpha: false });
  canvases.set(id, { ctx, imageData: ctx.createImageData(width, height) });
}

export function renderFrame(id, view) {
  const state = canvases.get(id);
  if (!state) {
    return;
  }
  // view is a MemoryView over WASM memory; slice() copies it out, which must
  // happen synchronously inside this call.
  state.imageData.data.set(view.slice());
  state.ctx.putImageData(state.imageData, 0, 0);
}

// --- vector-field arrow overlay (XY model) ---

const arrowCanvases = new Map();

export function initArrowCanvas(id, size) {
  const canvas = document.getElementById(id);
  if (!canvas) {
    return;
  }
  canvas.width = size;
  canvas.height = size;
  arrowCanvases.set(id, canvas.getContext('2d'));
}

export function clearArrows(id) {
  const ctx = arrowCanvases.get(id);
  if (ctx) {
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
  }
}

// angles is a MemoryView of gridSize*gridSize doubles (row-major); each is a
// representative spin direction for one cell of a coarse grid.
export function drawArrows(id, angles, gridSize) {
  const ctx = arrowCanvases.get(id);
  if (!ctx) {
    return;
  }

  // A MemoryView is not indexable; slice() copies it out to a Float64Array.
  const theta = angles.slice();
  const size = ctx.canvas.width;
  const cell = size / gridSize;
  const length = cell * 0.5;
  const head = cell * 0.2;

  ctx.clearRect(0, 0, size, size);
  ctx.strokeStyle = 'rgba(12, 12, 14, 0.55)';
  ctx.lineWidth = Math.max(1, cell * 0.05);
  ctx.lineCap = 'round';
  ctx.lineJoin = 'round';
  ctx.beginPath();

  for (let gy = 0; gy < gridSize; gy++) {
    for (let gx = 0; gx < gridSize; gx++) {
      const angle0 = theta[gy * gridSize + gx];
      const cx = (gx + 0.5) * cell;
      const cy = (gy + 0.5) * cell;
      const dx = Math.cos(angle0) * length;
      const dy = Math.sin(angle0) * length;
      const tipX = cx + dx / 2;
      const tipY = cy + dy / 2;
      const angle = Math.atan2(dy, dx);

      ctx.moveTo(cx - dx / 2, cy - dy / 2);
      ctx.lineTo(tipX, tipY);
      ctx.lineTo(tipX - head * Math.cos(angle - 0.5), tipY - head * Math.sin(angle - 0.5));
      ctx.moveTo(tipX, tipY);
      ctx.lineTo(tipX - head * Math.cos(angle + 0.5), tipY - head * Math.sin(angle + 0.5));
    }
  }

  ctx.stroke();
}

let running = false;

export function startLoop(dotNetRef) {
  if (running) {
    return;
  }
  running = true;
  const frame = async () => {
    if (!running) {
      return;
    }
    try {
      await dotNetRef.invokeMethodAsync('OnAnimationFrame');
    } catch {
      running = false;
      return;
    }
    if (running) {
      requestAnimationFrame(frame);
    }
  };
  requestAnimationFrame(frame);
}

export function stopLoop() {
  running = false;
}

function triggerDownload(url, filename) {
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
}

export function downloadCanvasPng(id, filename) {
  const canvas = document.getElementById(id);
  if (!canvas) {
    return;
  }
  canvas.toBlob((blob) => {
    const url = URL.createObjectURL(blob);
    triggerDownload(url, filename);
    URL.revokeObjectURL(url);
  }, 'image/png');
}

export function downloadText(filename, content) {
  const blob = new Blob([content], { type: 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  triggerDownload(url, filename);
  URL.revokeObjectURL(url);
}
