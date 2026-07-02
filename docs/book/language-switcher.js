# Shared language switcher for the oCAS mdBook sites.
#
# Each language book (en/, zh/) includes this via `additional-js` in its
# book.toml. It injects a "中文 / English" button into mdBook's top-right
# button bar that swaps the `/en/` and `/zh/` path segments while keeping
# the rest of the URL (version prefix, page, hash) intact.
#
# Deployed URL shapes this must handle:
#   /ocas/latest/en/introduction.html
#   /ocas/v0.10.0/zh/getting-started.html#rust
#   /en/404.html            (local `mdbook serve`)
;(function () {
  "use strict";

  function currentLang(path) {
    var m = path.match(/\/(en|zh)(\/|$)/);
    return m ? m[1] : null;
  }

  function swapLang(path, from, to) {
    // Replace only the first occurrence of /<from>/ to avoid touching a
    // version segment that happens to contain "en"/"zh".
    return path.replace("/" + from + "/", "/" + to + "/");
  }

  function setup() {
    var path = window.location.pathname;
    var cur = currentLang(path);
    if (!cur) return;

    var other = cur === "en" ? "zh" : "en";
    var label = other === "zh" ? "中文" : "English";
    var target = swapLang(path, cur, other);

    var anchor = document.createElement("a");
    anchor.href = target;
    anchor.className = "ocas-lang-switch";
    anchor.textContent = label;
    anchor.setAttribute("lang", other);
    anchor.title = "Switch language / 切换语言";

    // Prefer mdBook's right-button bar; fall back to the menu bar.
    var host =
      document.querySelector(".right-buttons") ||
      document.querySelector(".menu-bar") ||
      document.body;
    host.appendChild(anchor);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", setup);
  } else {
    setup();
  }
})();
