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
