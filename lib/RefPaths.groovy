class RefPaths {
    static String normalize(Object path) {
        if (path == null) {
            return null
        }
        def value = path.toString().trim()
        if (!value) {
            return null
        }
        while (value.endsWith('/') && !value.endsWith('://')) {
            value = value.substring(0, value.length() - 1)
        }
        return value
    }

    static String basename(Object path) {
        def value = normalize(path)
        if (!value) {
            return ''
        }
        def idx = Math.max(value.lastIndexOf('/'), value.lastIndexOf('\\'))
        return idx >= 0 ? value.substring(idx + 1) : value
    }

    static boolean isArchive(Object path) {
        def base = basename(path).toLowerCase()
        return base.endsWith('.tar.gz') || base.endsWith('.tgz') || base.endsWith('.tar')
    }
}
