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

// arrows is a MemoryView of count*3 doubles: (xPixel, yPixel, theta) per arrow.
// One arrow per lattice site, drawn only in the neighbourhood of vortices, so
// the per-site winding is resolved rather than averaged away.
export function drawVortexArrows(id, arrows, count, cellPixels) {
  const ctx = arrowCanvases.get(id);
  if (!ctx) {
    return;
  }

  // A MemoryView is not indexable; slice() copies it out to a Float64Array.
  const data = arrows.slice();
  const size = ctx.canvas.width;
  const length = cellPixels * 1.15;
  const head = cellPixels * 0.5;

  ctx.clearRect(0, 0, size, size);
  ctx.strokeStyle = 'rgba(10, 10, 12, 0.72)';
  ctx.lineWidth = Math.max(1, cellPixels * 0.11);
  ctx.lineCap = 'round';
  ctx.lineJoin = 'round';
  ctx.beginPath();

  for (let i = 0; i < count; i++) {
    const cx = data[i * 3];
    const cy = data[i * 3 + 1];
    const theta = data[i * 3 + 2];
    const dx = Math.cos(theta) * length;
    const dy = Math.sin(theta) * length;
    const tipX = cx + dx / 2;
    const tipY = cy + dy / 2;
    const angle = Math.atan2(dy, dx);

    ctx.moveTo(cx - dx / 2, cy - dy / 2);
    ctx.lineTo(tipX, tipY);
    ctx.lineTo(tipX - head * Math.cos(angle - 0.5), tipY - head * Math.sin(angle - 0.5));
    ctx.moveTo(tipX, tipY);
    ctx.lineTo(tipX - head * Math.cos(angle + 0.5), tipY - head * Math.sin(angle + 0.5));
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
