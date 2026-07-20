(function () {
  "use strict";

  const storagePrefix = "pymiesim-card-state:";

  function restoreCard(card) {
    if (!card.id || !card.matches("details.workflow-card")) {
      return;
    }

    const savedState = window.sessionStorage.getItem(storagePrefix + card.id);
    if (savedState === "open") {
      card.open = true;
    } else if (savedState === "closed") {
      card.open = false;
    }

    if (!card.dataset.statePersistenceBound) {
      card.addEventListener("toggle", function () {
        window.sessionStorage.setItem(storagePrefix + card.id, card.open ? "open" : "closed");
      });
      card.dataset.statePersistenceBound = "true";
    }
  }

  function scan(root) {
    if (root instanceof HTMLElement) {
      restoreCard(root);
      root.querySelectorAll("details.workflow-card").forEach(restoreCard);
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    scan(document.body);
    new MutationObserver(function (mutations) {
      mutations.forEach(function (mutation) {
        mutation.addedNodes.forEach(scan);
      });
    }).observe(document.body, { childList: true, subtree: true });
  });
})();
