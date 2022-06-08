def paramExists(name) {
  return params.containsKey(name) && params[name] != ""
}

def assertParamExists(name, description) {
  if (!paramExists(name)) {
    exit 1, "ERROR: Please provide a --${name} parameter ${description}"
  }
}

def getChild(parent, child) {
  if (child.contains("://") || java.nio.file.Paths.get(child).isAbsolute()) {
    child
  } else {
    parent.replaceAll('/[^/]*$', "/") + child
  }
}